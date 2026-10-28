"""
Convert GWTC event tables to asimov blueprints.
"""

from gwpy.table import EventTable
import astropy.units as u

import numpy as np
from gpstime import gpstime
import lalsimulation as lalsim
import lal
import conversions


table_path = "/data/wiay/GWTC-4.0/triggers/IGWN-GWTC4p0-0f954158d_720-SearchSummaryTable.hdf5"


blueprint_template = """
# IGR Catalogue Event Blueprint
kind: event
name: {name}
data:
  segment length: {seglen}
likelihood:
  psd length: {seglen}
  window length: {seglen}
  sample rate: {sample_rate}
interferometers:
{ifos}
priors:
  chirp mass:
    minimum: {chirp_mass[0]}
    maximum: {chirp_mass[1]}
event time: {event_time:.2f}
"""

def determine_event_name(event_time, prefix="GW"):
    """
    Determine the name of the event from an event time.
    """
    gt = gpstime.fromgps(event_time)
    time_string = gt.strftime('%y%m%d_%H%M%S')
    return prefix + time_string

def guess_chirp_mass_prior(chirp_mass, mass_ratio=1):
    """
    Guess a good chirp mass prior from a point estimate.
    """
    chirp_mass = int(chirp_mass)
    minimum = (chirp_mass - 15) if chirp_mass > 16 else 3
    maximum = chirp_mass + 15
    return minimum, maximum


def determine_segment_length(chirp_mass):
    """
    Determine a suitable segment length given the chirp mass.
    """
    chirp_mass *= lal.MSUN_SI
    m1, m2 = conversions.chirp_mass_q_to_components(chirp_mass, 1)
    length_1 = lalsim.SimIMRPhenomXASDuration(m1, m2,
                                            0.99, 0.99, 20)
    m1, m2 = conversions.chirp_mass_q_to_components(chirp_mass, 0.125)
    length_2 = lalsim.SimIMRPhenomXASDuration(m1, m2,
                                            0.99, 0.99, 20)
    length = int(length_1 if length_1 > length_2 else length_2)
    length += 1 # rollon
    length += 2 # post-trigger
    length += 0.03 #safety
    length = int(length)
    length = 1 << (length-1).bit_length()

    return length


def determine_sample_rate(chirp_mass):
    """
    Determine a sensible sampling rate.
    """
    M = 100.0  # will be ignored
    Ms = M * lal.MTSUN_SI
    m1 = M * 0.5
    m2 = M * 0.5
    complexQNM = lal.CreateCOMPLEX16Vector(1)
    lalsim.SimIMREOBGenerateQNMFreqV2(
        complexQNM,
        m1,
        m2,
        np.array([0.0, 0.0, 0.99]),
        np.array([0.0, 0.0, 0.99]),
        3, 3,
        1,
        lalsim.SEOBNRv4,
    )
    m1, m2 = conversions.chirp_mass_q_to_components(chirp_mass, 0.125)
    f = Ms * np.real(complexQNM.data[0]) / (lal.MTSUN_SI * (m1+m2) * 2.0 * np.pi)

    return 1 << (int(2*f)-1).bit_length()

def create_blueprint(chirp_mass, gpstime, ifos):
    name = determine_event_name(gpstime)
    chirp_prior = guess_chirp_mass_prior(chirp_mass)
    seglen = determine_segment_length(chirp_prior[0])
    sample_rate = max(determine_sample_rate(chirp_prior[0]),
                      determine_sample_rate(chirp_prior[1]))
    ifos_str = ""
    for ifo in ifos.split(","):
        ifos_str += f"  - {ifo}\n"
    return blueprint_template.format(name=name,
                                     chirp_mass=chirp_prior,
                                     event_time=gpstime,
                                     seglen=seglen,
                                     sample_rate=sample_rate,
                                     ifos=ifos_str,
                                     )



table = EventTable.read(table_path, path="search_summary")

far_threshold = (1/u.yr).to("Hz")
table_filtered = table.filter(f'p_Terrestrial < 0.5')
table_filtered = table_filtered.filter(f'far > {far_threshold.value}')

for data in table_filtered:

    name = determine_event_name(data['gps_time'], prefix="GW")
    with open(f"{name}.yaml", "w") as filename:
        filename.write(
            create_blueprint(data['mchirp'], data['gps_time'], ifos=data['instruments'])
        )

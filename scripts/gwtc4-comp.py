import os
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
from ledger import (
    load_ledger, get_events, approved_bilby_productions, open_pesummary_samples
)

if len(sys.argv) < 2:
    print(f"Usage: {sys.argv[0]} <event_name>", file=sys.stderr)
    sys.exit(1)

event = sys.argv[1]

web_path = os.environ.get("WEB_PATH", "/home/pe.o4/public_html/gwtc4-fogg")
comparison_web_path = os.environ.get("COMPARISON_WEB_PATH", "/home/pe.o4/public_html/gwtc4")
comparison_ledger = os.environ.get("COMPARISON_LEDGER", "/home/pe.o4/GWTC4/.asimov/ledger.yml")

waveforms = ["IMRPhenomXPHM-SpinTaylor",
             "SEOBNRv5PHM",
             "NRSur7dq4",
             "IMRPhenomXO4a"]

def reviewed_productions_for_event(event_name, ledger_path):
    data = load_ledger(ledger_path)
    events = get_events(data)
    if event_name not in events:
        return []
    return [p for p, _ in approved_bilby_productions(events[event_name])]


metafiles = [
    f"{web_path}/{event}/fogg-bilby-{waveform}/pesummary/samples/posterior_samples.h5"
    for waveform in waveforms
]

comparison_metafiles = [
    f"{comparison_web_path}/{event}/{production}/pesummary/samples/posterior_samples.h5"
    for production in reviewed_productions_for_event(event, comparison_ledger)
]

def tidy_plot(ax, xticks=False, yticks=False):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(left=False, bottom=False, labelsize=6)
    ax.grid(False)
    if not yticks:
        ax.set(yticklabels=[])
    if not xticks:
        ax.set(xticklabels=[])
    ax.tick_params(left=False, bottom=False, labelsize=6)
    ax.spines[['bottom', 'left', 'right', 'top']].set_visible(False)
    ax.grid(False)
    return ax

def plot_event_waveforms(
        name: str,
        metafiles: list,
        properties: list = [
            "chirp_mass",
            "luminosity_distance",
            "mass_1",
            "mass_2",
            "mass_ratio",
            "chi_eff",
            "chi_p",
            "comoving_distance",
        ]):


    fig = plt.figure(layout=None, figsize=(1.618*5, 0.33*len(properties)), dpi=300)
    gs = fig.add_gridspec(nrows=1+len(properties),
                          ncols=2+len(waveforms),
                          hspace=0.0, wspace=0.12)

    combined_ax = {}
    # Add the waveform names
    for column, waveform in enumerate(waveforms):
        ax = fig.add_subplot(gs[0, 2+column])
        ax.text(0,0.5, s=f"{waveform}", fontname="serif", verticalalignment="center", size=7)
        ax = tidy_plot(ax)
    
    # Add the property names
    for row, property in enumerate(properties):
        ax = fig.add_subplot(gs[1+row, 0])
        ax.text(0,0.5, s=f"{property}", fontname="serif", verticalalignment="center", size=7)
        ax = tidy_plot(ax)

        combined_ax[property] = fig.add_subplot(gs[1+row, 1])
        tidy_plot(combined_ax[property])
        
    for analysis_i, metafile in enumerate(metafiles):    
        if os.path.exists(metafile):
            posterior, prior = open_pesummary_samples(metafile)
            try:
                comp = [meta for meta in comparison_metafiles if waveforms[analysis_i] in meta][0]
                posterior_c, prior_c = open_pesummary_samples(comparison_metafiles[analysis_i])
            except IndexError:
                posterior_c = None

            for row, property in enumerate(properties):
                combined_ax[property].hist(posterior[property], histtype="step", bins=100)
                ax = fig.add_subplot(gs[1+row, 2+analysis_i])

                if posterior_c is not None:
                    ax.hist(posterior_c[property], histtype="step", bins=100, color="k", alpha=0.7)
                
                ax.hist(prior[property], histtype="stepfilled", color="k", alpha=0.3, bins=100)
                ax.hist(posterior[property], histtype="step", bins=100)

                ax = tidy_plot(ax)


    return fig

fig = plot_event_waveforms(event, metafiles)
output_path = os.path.join(web_path, event, "quickview.png")
os.makedirs(os.path.dirname(output_path), exist_ok=True)
fig.savefig(output_path)

"""Automated checks of gravitational-wave parameter estimation results."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("gw-catalog-tools")
except PackageNotFoundError:  # running from a source tree which isn't installed
    __version__ = "unknown"

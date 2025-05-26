from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("vmwhere")
except PackageNotFoundError:
    __version__ = "unknown"

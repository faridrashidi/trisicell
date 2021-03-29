import datetime
from platform import python_version

from trisicell import settings

_VERBOSITY_LEVELS_FROM_STRINGS = {
    "error": 0,
    "warn": 1,
    "info": 2,
    "hint": 3,
    "debug": 4,
}


def info(*args, **kwargs):
    return msg(*args, v="info", **kwargs)


def error(*args, **kwargs):
    args = ("Error:",) + args
    return msg(*args, v="error", **kwargs)


def warn(*args, **kwargs):
    args = ("WARNING:",) + args
    return msg(*args, v="warn", **kwargs)


def hint(*args, **kwargs):
    return msg(*args, v="hint", **kwargs)


def debug(*args, **kwargs):
    args = ("DEBUG:",) + args
    return msg(*args, v="debug", **kwargs)


def msg(*msg, v, time=False, end="\n"):
    if v is None:
        v = 4
    if isinstance(v, str):
        v = _VERBOSITY_LEVELS_FROM_STRINGS[v]
    if v == 3:
        msg = ("-->",) + msg
    if v >= 4:
        msg = ("   ",) + msg
    if _settings_verbosity_greater_or_equal_than(v):
        if len(msg) > 0:
            if time:
                msg = (f"[{datetime.datetime.now().strftime('%m/%d %H:%M:%S')}]",) + msg
            _write_log(*msg, end=end)


def _write_log(*msg, end="\n"):
    from trisicell.settings import logfile

    if logfile == "":
        print(*msg, end=end, flush=True)
    else:
        out = ""
        for s in msg:
            out += f"{s} "
        with open(logfile, "a") as f:
            f.write(out + end)
            f.flush()


def _settings_verbosity_greater_or_equal_than(v):
    from trisicell.settings import verbosity

    if isinstance(verbosity, str):
        settings_v = _VERBOSITY_LEVELS_FROM_STRINGS[verbosity]
    else:
        settings_v = verbosity
    return settings_v >= v


def _check_if_latest_version():
    from trisicell import __version__

    latest_version = timeout(
        get_latest_pypi_version, timeout_duration=2, default="0.0.0"
    )
    if __version__.rsplit(".dev")[0] < latest_version.rsplit(".dev")[0]:
        warn(
            "There is a newer trisicell version available on PyPI:\n",
            "Your version: \t\t",
            __version__,
            "\nLatest version: \t",
            latest_version,
        )


def _get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def print_version():
    from trisicell import __version__

    _write_log(
        f"Running trisicell {__version__} "
        f"(python {python_version()}) on {_get_date_string()}.",
    )
    # _check_if_latest_version()

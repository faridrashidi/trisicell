import datetime
from platform import python_version

import termcolor

_VERBOSITY_LEVELS_FROM_STRINGS = {
    "error": 0,
    "warn": 1,
    "info": 2,
    "hint": 3,
    "debug": 4,
}


def error(*args, **kwargs):
    """Write errors message."""
    args = ("Error:",) + args
    msg(*args, v="error", **kwargs)
    raise RuntimeError


def warn(*args, **kwargs):
    """Write warnings message."""
    args = ("WARNING:",) + args
    msg(*args, v="warn", **kwargs)


def info(*args, **kwargs):
    """Write info message."""
    msg(*args, v="info", **kwargs)


def hint(*args, **kwargs):
    """Write hints message."""
    msg(*args, v="hint", **kwargs)


def debug(*args, **kwargs):
    """Write debugs message."""
    args = ("DEBUG:",) + args
    msg(*args, v="debug", **kwargs)


def msg(*msg, v, time=False, color=None, end="\n"):
    r"""Write message to logging output.

    Parameters
    ----------
    v : :obj:`int`
        {'error', 'warn', 'info', 'hint', 'debug'} or int
        0/'error', 1/'warn', 2/'info', 3/'hint', 4/'debug', 5, 6...
    time : :obj:`bool`, optional
        Print date and time, by default False
    end : :obj:`str`, optional
        End character, by default "\n"
    """
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
            _write_log(*msg, color=color, end=end)


def _write_log(*msg, color=None, end="\n"):
    from trisicell.settings import logfile

    if logfile == "":
        if color is not None:
            termcolor.cprint(*msg, color, attrs=["bold"], end=end, flush=True)
        else:
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


def _get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def print_version():
    """Print version."""
    from trisicell import __version__

    _write_log(
        f"Running trisicell {__version__} "
        f"(python {python_version()}) on {_get_date_string()}.",
    )

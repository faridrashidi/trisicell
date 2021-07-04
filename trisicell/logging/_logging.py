import datetime
from platform import python_version

_VERBOSITY_LEVELS_FROM_STRINGS = {
    "error": 0,
    "warn": 1,
    "info": 2,
    "hint": 3,
    "debug": 4,
}


def info(*args, **kwargs):
    """[summary].

    Returns
    -------
    [type]
        [description]
    """
    msg(*args, v="info", **kwargs)


def error(*args, **kwargs):
    """[summary].

    Returns
    -------
    [type]
        [description]
    """
    args = ("Error:",) + args
    msg(*args, v="error", **kwargs)
    exit(1)


def warn(*args, **kwargs):
    """[summary].

    Returns
    -------
    [type]
        [description]
    """
    args = ("WARNING:",) + args
    msg(*args, v="warn", **kwargs)


def hint(*args, **kwargs):
    """[summary].

    Returns
    -------
    [type]
        [description]
    """
    msg(*args, v="hint", **kwargs)


def debug(*args, **kwargs):
    """[summary].

    Returns
    -------
    [type]
        [description]
    """
    args = ("DEBUG:",) + args
    msg(*args, v="debug", **kwargs)


def msg(*msg, v, time=False, end="\n"):
    r"""[summary].

    Parameters
    ----------
    v : [type]
        [description]
    time : bool, optional
        [description], by default False
    end : str, optional
        [description], by default "\n"
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


def _get_date_string():
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M")


def print_version():
    """[summary]."""
    from metools import __version__

    _write_log(
        f"Running metools {__version__} "
        f"(python {python_version()}) on {_get_date_string()}.",
    )

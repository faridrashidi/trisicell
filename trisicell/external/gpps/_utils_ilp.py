__author__ = "Simone Ciccolella"
__date__ = "11/30/21"


def expand_name(s, max_gains, max_losses):
    positive_names = [s + "+" + str(i) for i in range(max_gains)]
    negative_names = [s + "-" + str(i) for i in range(max_losses)]
    return positive_names + negative_names

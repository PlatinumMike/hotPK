# -*- coding: utf-8 -*-
"""
Compute S functions.
Requires the mpmath package because it has the meijerg, but also because it
needs the arbitrary precision at large arguments.
"""

import mpmath as mp

mp.mp.dps = 150  # set precision


def getS1(eps: int, x: float):
    return 0.5j * mp.meijerg([[], []], [[0, 0], [1 / 2]], x * x) + eps * 0.5 * (
        4 * abs(x) * mp.hyper([], [3 / 2, 3 / 2], x * x) - mp.sqrt(mp.pi) * mp.hyper([], [1 / 2, 1], x * x)
    )


def getS2(eps: int, x: float):
    return eps * x * mp.meijerg([[], []], [[0, 0], [-1 / 2]], x * x) + mp.sign(x) * 1.0j * (
        2 * mp.sqrt(mp.pi) * abs(x) * mp.hyper([], [1, 3 / 2], x * x) - mp.hyper([], [1 / 2, 1 / 2], x * x)
    )


def getS3(eps: int, x: float):
    return -1.0j * mp.meijerg([[], []], [[0, 1], [1 / 2]], x * x) + 2 * eps * (
        abs(x) * mp.hyper([], [1 / 2, 3 / 2], x * x) - mp.sqrt(mp.pi) * x * x * mp.hyper([], [3 / 2, 2], x * x)
    )

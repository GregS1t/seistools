#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Useful time functions developped by ETH




"""
import datetime as dt
from obspy import UTCDateTime
import numpy as np

# Earth and Mars days in seconds
SECONDS_PER_EARTH_DAY = 86400.

# Sol0 start - Insight landing.
# Value is from the JPL URL below for Sol 0:
# https://naif.jpl.nasa.gov/cgi-bin/chronos_nsyt.pl?setup=nsyttime
REFERENCE_TIME_INSIGHT_LANDING = UTCDateTime("2018-11-26T05:10:50.336037Z")

# Sol-001 and Sol-002 start times to compute one Martian day in seconds.
# Cannot use landing time because Sol-000 lasted shorter.
SOL01_START_TIME = UTCDateTime("2018-11-27T05:50:25.580014Z")
SOL02_START_TIME = UTCDateTime("2018-11-28T06:30:00.823990Z")

# Compute one martian day in seconds with a microsecond correction
# to avoid falling into previous or next sol instead of current
# one at sol midnights
SECONDS_PER_MARS_DAY = SOL02_START_TIME - SOL01_START_TIME - 0.000005


def lmst2utc(lmst_time, sol0_start_utc=REFERENCE_TIME_INSIGHT_LANDING):
    """
    Returns corresponding UTC date/time for a given LMST value. LMST can
    be a float or UTCDateTime instance. LMST variable should take Linux
    epoch as basis.

    Return value is an UTCDateTime instance.
    """
    _mars_to_earth = float(lmst_time) / SECONDS_PER_EARTH_DAY + 1

    _utc_time = UTCDateTime(
        _mars_to_earth * SECONDS_PER_MARS_DAY + float(sol0_start_utc))

    return _utc_time


def utc2lmst(utc_time, sol0_start_utc=REFERENCE_TIME_INSIGHT_LANDING,
             sol_dtype='int'):
    """
    Convert UTC to LMST. Default sol-0 time is InSight landing time
    in UTC. Returned LMST counts from Linux epoch; date value showing
    the sol number.

    Return value is a tuple of LMST as UTCDateTime instance and sol number.
    Sol number can be integer or float. If float, it includes decimal
    fractions of sol as well.
    """
    # Cast to UTCDateTime, if datetime is given. This is useful to avoid
    # long type casting statements while plotting
    if isinstance(utc_time, dt.datetime):
        utc_time = UTCDateTime(utc_time)

    _elapsed_mars_days = (utc_time - sol0_start_utc) / SECONDS_PER_MARS_DAY

    _time = UTCDateTime((_elapsed_mars_days - 1) * SECONDS_PER_EARTH_DAY)

    # Return a tuple with local Mars time as UTCDateTime and sol number
    if sol_dtype == 'float' or sol_dtype is None:
        return _time, _elapsed_mars_days
    else:
        return _time, np.int(np.floor(_elapsed_mars_days))


def utc2sol(utc_time, sol0_start_utc=REFERENCE_TIME_INSIGHT_LANDING,
            sol_dtype='int'):
    """
    A short hand call to utc2lmst method to get sol number. Return value
    is the sol number as integer or float, depending on sol_dtype.

    sol_dtype is type of string.
    """
    _, sol_number = utc2lmst(
        utc_time=utc_time, sol0_start_utc=sol0_start_utc, sol_dtype=sol_dtype)

    return sol_number


def sol_span_in_utc(sol, sol0_start_utc=REFERENCE_TIME_INSIGHT_LANDING):
    """
    Returns start and end times in UTC for a given sol.
    """
    utc_representation = \
        UTCDateTime(sol * SECONDS_PER_MARS_DAY) + float(sol0_start_utc)

    return utc_representation, utc_representation + SECONDS_PER_MARS_DAY


def parse_location_channel_codes(channel):
    """
    Separate location and channel codes in a string that is given in
    a form of '58.BZC'
    """
    try:
        dot_index = channel.index('.')
        return channel[0:dot_index], channel[dot_index+1:]

    except ValueError:
        raise ValueError('Channel name should include a dot: %s' % channel)

def dayplot_set_x_ticks(starttime, endtime, ax=None, print_sol=False,
                        show_last_tick=True, plot_doy=False,
                        sol_time_format=None, **kwargs):
    """
    Sets time label xticks for the Sol plot. If print_sol is False, labels
    are generated for UTC.
    """
    hour_ticks = []
    tick_labels = []
    hour_ticks_minor = []

    interval = np.abs(endtime - starttime)
    interval_h = interval / 3600.

    # Ticks start at the first hour of waveform time span
    tick_start = UTCDateTime(
        starttime.year, starttime.month, starttime.day, starttime.hour)

    # Tick label in every 3 marks. Enough for sol spectrograms
    step = 3 * 3600
    step_h = step / 3600.

    if print_sol:
        # Start from 1 to start from midnight in LMST
        start_index = 1

        for _ihour in np.arange(start_index, interval_h, step_h):
            hour_tick = tick_start + _ihour * 3600.
            hour_ticks.append(float(hour_tick))

            if sol_time_format is None:
                # Get sol number to construct the tick label
                _lmst_start, _sol = utc2lmst(lmst2utc(hour_tick))
                _tick_label = 'Sol-' + str(_sol).zfill(3) + '\n' + \
                              UTCDateTime(hour_tick).strftime('%H:%M')
            else:
                _tick_label = UTCDateTime(hour_tick).\
                    strftime(sol_time_format)

            tick_labels.append(_tick_label)

        if show_last_tick:
            hour_tick = tick_start + (_ihour + step_h) * 3600.
            hour_ticks.append(float(hour_tick))

            if sol_time_format is None:
                # Get sol number to construct the tick label
                _lmst_start, _sol = utc2lmst(lmst2utc(hour_tick))
                _tick_label = 'Sol-' + str(_sol).zfill(3) + '\n' + \
                              UTCDateTime(hour_tick).strftime('%H:%M')
            else:
                _tick_label = UTCDateTime(hour_tick). \
                    strftime(sol_time_format)

            tick_labels.append(_tick_label)

    else:
        start_index = 1

        for _ihour in np.arange(start_index, interval_h, step_h):
            hour_tick = tick_start + _ihour * 3600.
            hour_ticks.append(float(hour_tick))

            # Always plot day of year to first label
            if _ihour == start_index:
                if plot_doy:
                    tick_labels.append(
                        UTCDateTime(hour_tick).strftime(
                            '%H:%M%n%Y-%m-%d%n[DOY %j]'))
                else:
                    tick_labels.append(
                        UTCDateTime(hour_tick).strftime(
                            '%Y-%m-%d%n%H:%M'))

            # Otherwise, plot day of year when day changes
            else:
                _previous_tick_time = UTCDateTime(
                    tick_start + (_ihour - step_h) * 3600.)

                if UTCDateTime(hour_tick).julday != \
                        _previous_tick_time.julday:
                    if plot_doy:
                        tick_labels.append(
                            UTCDateTime(hour_tick).strftime(
                                '%H:%M%n%Y-%m-%d%n[DOY %j]'))
                    else:
                        tick_labels.append(
                            UTCDateTime(hour_tick).strftime(
                                '%H:%M'))
                else:
                    tick_labels.append(
                        UTCDateTime(hour_tick).strftime(
                            '%H:%M'))

    # Minor ticks - no labels
    for _ihour in np.arange(0, interval_h + 1, 1):
        hour_tick = tick_start + _ihour * 3600.
        hour_ticks_minor.append(float(hour_tick))

    # If axes is given, set the ticks and labels.
    if ax:
        ax.set_xticks(hour_ticks)
        ax.set_xticks(hour_ticks_minor, minor=True)
        ax.set_xticklabels(tick_labels, **kwargs)

        ax.set_xlim(float(starttime), float(endtime))

    # Return the tick locations and major tick labels
    return {'major-ticks': hour_ticks, 'minor-ticks': hour_ticks_minor,
            'labels': tick_labels}



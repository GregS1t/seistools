#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 13 14:26:23 2019

@author  : Grégory Sainton
@email   : sainton`at`ipgp.fr
@purpose : Functions to handle julian dates

CNES if often using Julian date format.
Here are to simple function to convert `normal` date format
    Julian date -> date
    date -> Julian date

Julian date format: YYYY-DDD
Normal date format: YYYY-MM-DD


"""

import calendar
from datetime import date

def from_juliandate(year, juldate):
    """
    Code to convert Julian date to classical MMDDYYYY
    date format

    Parameters
    ----------
    y : int
        DESCRIPTION.
    jd : int
        DESCRIPTION.

    Returns
    -------
        date - string - YYYY-MM-DD format

    """
    month = 1
    while juldate - calendar.monthrange(year, month)[1] > 0 and month <= 12:
        juldate = juldate - calendar.monthrange(year, month)[1]
        month = month + 1

    return f"{year}-{month}-{juldate}"


def to_juliandate(inputdate):
    """
    Function to convert date (YYYY-MM-DD) into a Julian format date
    (YYYYDDD)

    Parameters
    ----------
    date : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    splitteddate = inputdate.split("-")
    if len(splitteddate)==3:
        newdate = date(int(splitteddate[0]), int(splitteddate[1]),
                    int(splitteddate[2]))
        date1 = date(newdate.year, 1, 1)
        jday = newdate - date1 
        juliandate = f"{newdate.year}-{jday.days+1:03d}"
        return(juliandate)
    else:
        raise Exception("Sorry, your date must be with the following format: YYYY-MM-DD")


def main():
    """
    Main program: Example of execution.
    """
    print("2 examples to use these functions: \n")

    # from_juliandate(year, juldate)
    new_date = from_juliandate(2022, 150)
    print(new_date)

    #to_juliandate(inputdate)
    new_julian_date = to_juliandate("2022-01-13")
    print(new_julian_date)

if __name__ == "__main__":
    main()
#!/usr/bin/env python
# -*- coding: utf-8 -*-


class ENAFetchFail(Exception):
    """Raised when max number of attempts to fetch data from the API is reached"""
    pass


class ENAFetch401(Exception):
    """Raised on authentication error"""
    pass


class ENAFetch204(Exception):
    """Raised when no data found for the given accession"""
    pass


class NoDataError(Exception):
    """Raised when no run, assembly or study data is found"""
    pass

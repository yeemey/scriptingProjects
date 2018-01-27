#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 15:51:17 2018

functions for mortgage-related calculations

@author: ymseah
"""

def pmi_rate(monthly_pmi, initial_loan_amt):
    pmi_rate = (monthly_pmi * 12)/initial_loan_amt
    return pmi_rate

def monthly_pmi(pmi_rate, initial_loan_amt):
    monthly_pmi = (pmi_rate * initial_loan_amt)/12
    return monthly_pmi

def monthly_principal_interest(monthly_payment, ann_interest_as_percent, initial_loan_amt, which_month=1):
    '''
    monthly_payment = monthly payment of principal + interest only, not including taxes, homeowners insurance, PMI, etc.
    '''
    mth_counter = 1
    monthly_interest_rate = (ann_interest_as_percent/100)/12
    remaining_principal = initial_loan_amt
    mth_counter = 1
    while mth_counter <= which_month:
        month_interest = monthly_interest_rate * remaining_principal
        month_principal = monthly_payment - month_interest
        remaining_principal -= month_principal
        mth_counter += 1
    return month_principal, month_interest

def principal_interest_to_date(monthly_payment, ann_interest_as_percent, initial_loan_amt, which_month=1):
    '''
    monthly_payment = monthly payment of principal + interest only, not including taxes, homeowners insurance, PMI, etc.
    '''
    principal = 0
    interest = 0
    mth_counter = 1
    while mth_counter <= which_month:
        month_principal, month_interest = monthly_principal_interest(monthly_payment, ann_interest_as_percent, initial_loan_amt, which_month)
        principal += month_principal
        interest += month_interest
        mth_counter += 1
    return principal, interest

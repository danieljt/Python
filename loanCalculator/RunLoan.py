from Loancalculator import *

# Set values
capital = 2000000
interest = 3.0
years = 25

# calculate annuity loan
loan = Series(capital,interest,years)
loan.tabulate()

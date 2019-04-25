from scitools.std import *

class Loan:
    """
    Superclass for defining a loan. Program sets the interest, the
    amount of years for the loanto run, and the initial loan
    amount. Currently, the loancalculator does not take into account
    monthly and startup fees.  

    -------------------------------------------------
    USAGE:
    loan = Loan(capital, interest, years)
    -------------------------------------------------
    """
    def __init__(self, capital, interest, years):
        self.capital = capital
        self.interest = interest
        self.years = years
        self.months = years*12

        p_factor_y = interest/100.
        g_factor_y = 1 + p_factor_y
        g_factor_m = g_factor_y**(1./12)
        p_factor_m = g_factor_m - 1

        self.g_factor_m = g_factor_m
        self.p_factor_m = p_factor_m

class Annuity(Loan):
    """
    Class for running an annuity loan for the tabulated values
    in the loan. the class takes the parameters from the
    superclass Loan and calculates the loan.

    ---------------------------------------------------
    USAGE:
    loan = Annuity(capital, interest, years)
    -------------------------------------------------
    """
    def tabulate(self):
        """
        Creates a table of the monthly payback plan
        for the annuity loan.

        ---------------------------------------------
        USAGE:
        loan = Annuity(capital, interest, years)
        loan.tabulate()
        ---------------------------------------------
        """

        # Sets variables for the loan from the superclass
        capital = self.capital
        months = self.months
        percent = self.p_factor_m
        growth = self.g_factor_m

        # Calculates the installments
        installment = capital*growth**months*(growth - 1)/(growth**months - 1)

        # Set starting values
        loan_left = capital
        total_rate_paid = 0
        total_cost = 0 
        i = 0

        # Prints the table
        print 20*'%%%'
        print '|','CALCULATED MONTHLY PAYMENT PLAN BY AN ANNUITY LOAN'
        print 20*'%%%'
        print '|','%12s %12s %12s %12s' % \
            ('Base amount','Monthly','Interest','Loan paid'),'|'
        while i <= months - 1:
            rate_paid = percent*loan_left
            loan_paid = installment - rate_paid
            total_rate_paid += rate_paid
            total_cost += installment

            print '|','%12.2f %12.2f %12.2f %12.2f' %\
                (loan_left, installment, rate_paid, loan_paid),'|'
            loan_left -= loan_paid
            i += 1
        print 20*'---'
        print '|','TOTAL AMOUNT: ', total_cost
        print '|','TOTAL COST  : ', total_rate_paid
        print 20*'---'

class Series(Loan):
    """
    Class for running a Series loan for the tabulated values
    in the loan. the class takes the parameters from the
    superclass Loan and calculates the loan.

    ---------------------------------------------------
    USAGE:
    loan = Annuity(capital, interest, years)
    -------------------------------------------------
    """
    def tabulate(self):
        """
        Creates a table of the monthly payback plan
        for the series loan.

        ---------------------------------------------
        USAGE:
        loan = Series(capital, interest, years)
        loan.tabulate()
        ---------------------------------------------
        """

        # Sets the loan variables from the superclass
        capital = self.capital
        months = self.months
        percent = self.p_factor_m
        growth = self.g_factor_m

        # Starting values
        loan_left = capital
        monthly = float(capital)/months
        total_rate_paid = 0
        total_cost = 0 
        i = 0

        #Prints the payback plan
        print 20*'%%%'
        print '|','CALCULATED MONTHLY PAYMENT PLAN BY A SERIES LOAN'
        print 20*'%%%'
        print '|','%12s %12s %12s %12s' % \
            ('Base amount', 'Monthly', 'Interest', 'Loan Paid'),'|'
        while i <= months - 1:
            rate_paid = percent*loan_left
            installment = monthly + rate_paid
            total_cost += installment
            total_rate_paid += rate_paid

            print '|','%12.2f %12.2f %12.2f %12.2f' %\
                (loan_left, installment, rate_paid, monthly),'|'
            loan_left -= installment
            i += 1
        print 20*'---'
        print '|','TOTAL AMOUNT :', total_cost
        print '|','TOTAL COST   :', total_rate_paid
        print 20*'---'

    
if __name__=='__main__':
    """
    Test program for running the annuity and series loan calculators
    """
    capital = 2000000
    interest = 3.0
    years = 3
    
    loan = Annuity(capital, interest, years)
    loan.tabulate()

    loan = Series(capital, interest, years)
    loan.tabulate()

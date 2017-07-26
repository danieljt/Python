def iterator(filename):
    """
    Iterates the numbers 0 - 9 and solves the equation
    _ _ _ * _ _ - _ _ _ * _ _ =
    where the numbers are added over the _ signs. 
    """
    f = open(filename, 'w') 
    counter = 0;
    remaining = [0,1,2,3,4,5,6,7,8,9]
    for n0 in remaining:
        current0 = list(remaining)
        current0.remove(n0)
        for n1 in current0:
            current1 = list(current0)
            current1.remove(n1)
            for n2 in current1:
                current2 = list(current1)
                current2.remove(n2)
                for n3 in current2:
                    current3 = list(current2)
                    current3.remove(n3)
                    for n4 in current3:
                        current4 = list(current3)
                        current4.remove(n4)
                        for n5 in current4:
                            current5 = list(current4)
                            current5.remove(n5)
                            for n6 in current5:
                                current6 = list(current5)
                                current6.remove(n6)
                                for n7 in current6:
                                    current7 = list(current6)
                                    current7.remove(n7)
                                    for n8 in current7:
                                        current8 = list(current7)
                                        current8.remove(n8)
                                        for n9 in current8:
                                            res = subtract(n0, n1, n2,n3, n4, \
                                                          n5, n6, n7, n8, n9)

                                            f.write("%3s" % counter)
                                            f.write("%5s" % res[0])
                                            f.write("%4s" % res[1])
                                            f.write("%5s" % res[2])
                                            f.write("%3s" % res[3])
                                            f.write("%7s" % res[4])
                                            f.write("\n")

                                            counter += 1
    f.close()
                                            
def subtract(n0, n1, n2, n3, n4, n5, n6, n7, n8, n9):
    """
    Create the numbers, subtract and return everything
    """
    num1 = int(str(n0) + str(n1) + str(n2))
    num2 = int(str(n3) + str(n4))
    num3 = int(str(n5) + str(n6) + str(n7))
    num4 = int(str(n8) + str(n9))

    return num1, num2, num3, num4, num1*num2 - num3*num4

def findmin(filename):
    f = open(filename, 'r')

    # Lists for min values
    minlist = list()
    counter = 0

    for line in f:
        digits = line.split()
        minimum = int(digits[5])
        if(minimum == 0):
            minlist.append(digits)
            counter += 1
            
    f.close()
    print 20*'--'
    print "RESULTS:"
    print "There are %d different combinations that work" % counter
    print ""
    for i in range(len(minlist)):
        print "%3s %10s %4s %3s %4s %3s %2s" % (i, minlist[i][0], minlist[i][1],
                                               minlist[i][2], minlist[i][3],
                                               minlist[i][4], minlist[i][5])
    print 20*'--'
    
def main():
    filename = 'workfile'
    iterator(filename)
    findmin(filename)

if __name__ == '__main__':
    main()

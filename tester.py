import csv
import os.path
import sys



def tester1():
    X = []
    Y = []
    cl = []
    Cp = []

    relpath = sys.path[0] + '/tests/test_data/NACA0012.csv'
    with open(relpath,'r') as csvfile:
        spamreader = csv.reader(csvfile,delimiter=',')
        testData = [None]*4
        counter = 0
        for row in spamreader:
            if (X == []):
                X = row
            elif (Y == []):
                Y = row
            elif (cl == []):
                cl = row
            elif (Cp == []):
                Cp = row

        print(X)

def main():
    tester1()


if __name__ == "__main__":
    main()

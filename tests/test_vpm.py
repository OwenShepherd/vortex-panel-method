import os.path
import sys

def func(x):
    return (x+1)


def test_answer():
    assert func(3) == 5


def test_NACA0012():
    script_path = sys.path[0]
    relpath = '/tests/test_data/NACA0012.csv'
    abspath = script_path + relpath

    with open(abspath,'r') as csvfile:
        spamreader = csv.reader(csvfile,delimiter=',')
        







class TestClass(object):
    def test_one(self):
        x = 'this'
        assert 'h' in x

    def test_two(self):
        x = 'hello'
        assert hasattr(x, 'check')

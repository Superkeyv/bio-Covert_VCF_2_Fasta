from unittest import TestCase
import gzip
import time

class Test(TestCase):
    def test_line_process(self):
        self.fail()



def get_file_line_num(f):
    '''
    用于确认文件有多少行
    :param f: 所打开的文件名
    :return: 该文件存在的行数
    '''
    line_num = 0

    start = time.time_ns()
    while (True):
        buf = f.read(4 * 1024 * 1024)
        if (buf == ''):
            break

        line_num += buf.count('\n')

    time_used = time.time_ns() - start

    return line_num, time_used


class Test(TestCase):
    def test_get_file_line_num(self):
        f1 = open('a.txt', 'rt')
        f2 = gzip.open('a.txt.gz', 'rt')

        l1, t1 = get_file_line_num(f1)
        l2, t2 = get_file_line_num(f2)

        print("l1:{},t1:{}".format(l1, t1))
        print("l2:{},t2:{}".format(l2, t2))

        f1.close()
        f2.close()



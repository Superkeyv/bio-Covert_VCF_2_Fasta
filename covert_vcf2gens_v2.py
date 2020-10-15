# 这是并行化的数据处理实现，joblib

import argparse
import pandas as pd
import os
from progressbar import ProgressBar
from ParallelLineProcess.LinePrcessor import ParallelLine
import multiprocessing as mp

print_prefix = "STATUS -- "

# 碱基配对字典
BasePair = {'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C'}


def get_file_line_num(fname):
    '''
    用于确认文件有多少行
    :param f: 所打开的文件名
    :return: 该文件存在的行数
    '''
    line_num = 0
    with open(fname, 'r') as f:
        while (True):
            buf = f.read(4 * 1024 * 1024)
            if (buf == ''):
                break

            line_num += buf.count('\n')

        f.close()
    return line_num


def proc_col_of_pdchunk(data):
    '''
    用于处理chunk的一个子任务
    :param data: 传入的参数。元组类型(col_name, ck_col)
    :param ck_col: 需要处理的一个chunk列
    :param col_name: 该列对应的名称
    :return:
    '''
    col_name = data[0]
    ck_col = data[1]

    # 追加写文件
    f = open('tmp/{}.genseries'.format(col_name), 'a+')
    for GT in ck_col.values:
        f.write(GT)
    f.close()  # 文件变动写入磁盘


def covert_csv2phylib(opt):
    '''
    根据转化得到的基因型拼接得到phylib形式的基因序列
    文件以行为单位，每一行代表一个样本，比如s1
    :param opt: 使用的参数
    :return:
    '''
    with open(opt.file + '.gen_series.phy', 'w') as genseries_file:
        print(print_prefix + '正在进行基因序列生成')

        opt.csv_size = os.path.getsize(opt.file + '.csv')
        print(print_prefix + "csv文件大小: {} bytes".format(opt.csv_size))

        # 获取基因长度信息
        csv_reader = pd.read_csv(opt.file + '.csv', iterator=True, chunksize=2)
        title = csv_reader.read(1)
        opt.sample_size = len(title.columns) - 2
        csv_reader.close()
        opt.gen_length = get_file_line_num(opt.file + '.csv')

        if (opt.genfile_type == 'phylib'):
            # 输出基因相关的数据
            genseries_file.write('{} {}\n'.format(opt.sample_size, opt.gen_length * 2))

        if (opt.with_limited_mem):
            csv_reader = pd.read_csv(opt.file + '.csv', iterator=True, chunksize=2)
            print(print_prefix + '样本总数为:{}, 基因长度:{}'.format(opt.sample_size, opt.gen_length))

            # 样本列表
            sample_list = []
            for i in range(1, opt.sample_size + 1):
                sample_list.append('s{}'.format(i))

            # 同时创建对应文件
            if (not os.path.exists('tmp')):
                os.mkdir('tmp')

            # 考虑到内存太小，这里使用缓存文件作为中间处理
            for colname in sample_list:
                f = open('tmp/{}.genseries'.format(colname), 'w')
                f.close()

            # 处理进度相关
            progressbar = ProgressBar(maxval=opt.csv_size)
            progressbar.start()
            processed_size = 0

            # 启动进程池
            pool = mp.Pool(opt.parallel_jobs)

            # 按块读取csv文件，同时输出到对应基因序列中
            csv_reader = pd.read_csv(opt.file + '.csv', iterator=True, chunksize=opt.with_csv_chunksize)
            for ck in csv_reader:
                # 估计当前的进度
                processed_size += ck.size * 3
                progressbar.update(int(processed_size))

                if (opt.parallel_jobs > 1):

                    # 将数据包装上列明
                    arg_list = []
                    for sample_name in sample_list:
                        arg_list.append((sample_name, ck[sample_name]))

                    pool.imap(proc_col_of_pdchunk, arg_list)

                    # Parallel(n_jobs=opt.parallel_jobs)(
                    #     delayed(proc_col_of_pdchunk)(ck[sample_name], sample_name) for sample_name in sample_list
                    # )

                else:
                    # 将csv的块加载进内存，降低内存使用量
                    for sample_name in sample_list:
                        # 处理1个chunk的内容，将基因序列写入磁盘
                        ck_col = ck[sample_name]
                        # 追加写文件
                        f = open('tmp/{}.genseries'.format(sample_name), 'a+')
                        for GT in ck_col.values:
                            f.write(GT)
                        f.close()  # 文件变动写入磁盘

            # 处理完成
            progressbar.finish()

            # 等待各进程完成工作
            pool.close()
            pool.join()

            # 将所有基因文件合并，并清理tmp文件夹
            print(print_prefix + "开始合并文件")
            progressbar = ProgressBar(maxval=len(sample_list))
            progressbar.start()
            for i, colname in enumerate(sample_list):
                progressbar.update(i)
                sX_f = open('tmp/{}.genseries'.format(colname), 'r')

                if (opt.genfile_type == 'phylib'):
                    genseries_file.write('seq{}\t'.format(i + 1))
                genseries_file.write('{}\n'.format(sX_f.readline()))

            progressbar.finish()
            print(print_prefix + "文件合并完成")

        else:
            csv_reader = pd.read_csv(opt.file + '.csv')

            if (not hasattr(opt, 'sample_size')):
                opt.sample_size = len(csv_reader.columns) - 2
            print(print_prefix + '样本总数为:{}'.format(opt.sample_size))
            progressbar = ProgressBar(maxval=opt.sample_size)
            progressbar.start()

            for i, col_name in enumerate(csv_reader.columns[2:]):
                # 更新文件处理进度
                progressbar.update(i)

                col = csv_reader[col_name]

                for GT in col.values:
                    genseries_file.write(GT)
                genseries_file.write('\n')

            progressbar.finish()

        print(print_prefix + '基因序列生成完毕，序列总数:{}'.format(opt.sample_size))


ignore_ALTs = 0
loss_replace = '-'


def line_process(line):
    # 处理文件头
    if (line[0:2] == '##'):
        return

    if (line[0:2] == '#C'):
        # 这里是第一行，可以作为头部
        head = line.split()
        title = 'CHROM,POS'

        s_index = 0
        for i, col_name in enumerate(head):
            if (col_name[0] == 's'):
                s_index = i
                break

        for col in head[s_index:]:
            title += ',{}'.format(col)
        return title

    # 解析处理单行数据
    s_data = line.split()

    # 处理ALT特殊情况，当变异类型超过opt.ignore_ALTs时，忽略当前行的记录
    if (ignore_ALTs > 0 and ignore_ALTs < len(s_data[4].split(','))):
        return

    # 记录染色体，SNP位点位置
    row = s_data[0] + ',' + s_data[1]

    # 准备参考碱基与对应碱基
    REF = s_data[3]
    RREF = BasePair[REF]

    # 从9列到最后，进行遍历处理
    for i in range(9, len(s_data)):
        s = s_data[i][0:3]

        GT = ''  # 当前样本该位置的基因型
        if (s == './.'):
            GT = loss_replace * 2
        else:
            # 在没有缺失的情况下，解析基因型
            if (s[0] == '0'):
                GT += REF
            else:
                GT += RREF

            if (s[2] == '0'):
                GT += REF
            else:
                GT += RREF

        row += ',' + GT

    return row


def covert_vcf2csv(opt):
    '''
    根据参数设置，解析数据
    :param opt: 使用的参数
    :return:
    '''
    src = open(opt.file, 'r')
    out = open(opt.file + '.csv', 'w')

    pline = ParallelLine(n_jobs=opt.parallel_jobs, chunk_size=opt.with_csv_chunksize, show_process_status=True)
    pline.run_row(input_file=src, output_file=out, order=True, row_func=line_process, use_CRLF=True)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="本脚本用于处理vcf数据，最后生成基因序列",
                                     description='示例：python covert_vcf2genseries.py --file sample.vcf gen_series True '
                                                 '--ignore_ALTs 1 --with_csv_chunksize 10000 --parallel_jobs 4 ')
    parser.add_argument('--file', type=str, default='sample.vcf', help="需要处理的文件")
    parser.add_argument('--ignore_ALTs', type=int, default='0', help="变异类型超过该值则忽略该行。0表示不忽略")
    parser.add_argument('--loss_replace', type=str, default='-', help="基因型缺失的替代符号，需要用引号包括")
    parser.add_argument('--covert_vcf2csv', type=str2bool, default=True,
                        help="是否执行covert_vcf2csv部分。如果已经转化出csv文件，设置true可以节省大量时间。可选参数：True,False")
    parser.add_argument('--gen_series', type=str2bool, default=True,
                        help="是否生成基因序列，建议设置基因型缺失时的替代符号，默认用空格代替。可选参数：True,False")
    parser.add_argument('--with_limited_mem', type=str2bool, default=True,
                        help="内存的使用策略。True尽量少的利用内存(速度相对慢一些)，False尽量多的利用内存(速度相对快一些)。可选参数：True,False")
    parser.add_argument('--with_csv_chunksize', type=int, default=1000,
                        help="设置chunksize，控制csv分块读取时的行数。面对超大文件时可以有效降低内存使用。设置过大会延长加载时间")
    parser.add_argument('--parallel_jobs', type=int, default=4,
                        help="并行操作使用的进程数，设置1代表不并行。由于磁盘读取速度的限制。建议这里不要设置超过8")
    parser.add_argument('--genfile_type', type=str, default='phylib',
                        help="生成的基因序列文件格式。可选有'txt'、'phylib'")

    opt = parser.parse_args()

    ignore_ALTs = opt.ignore_ALTs
    loss_replace = opt.loss_replace

    print(opt)
    if (opt.covert_vcf2csv):
        covert_vcf2csv(opt)

    if (opt.gen_series):
        covert_csv2phylib(opt)
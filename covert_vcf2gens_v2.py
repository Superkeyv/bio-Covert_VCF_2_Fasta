# 这是并行化的数据处理实现，joblib

import argparse
import pandas as pd
import os
from progressbar import ProgressBar
from ParallelLineProcess.LinePrcessor import ParallelLine
import multiprocessing as mp
import shutil

print_prefix = "STATUS -- "
tmp_file_dir = 'tmp'

# IUPAC碱基编码表
amb = {"*": "-", "A": "A", "C": "C", "G": "G", "N": "N", "T": "T",
       "*A": "a", "*C": "c", "*G": "g", "*N": "n", "*T": "t",
       "AC": "M", "AG": "R", "AN": "a", "AT": "W", "CG": "S",
       "CN": "c", "CT": "Y", "GN": "g", "GT": "K", "NT": "t",
       "*AC": "m", "*AG": "r", "*AN": "a", "*AT": "w", "*CG": "s",
       "*CN": "c", "*CT": "y", "*GN": "g", "*GT": "k", "*NT": "t",
       "ACG": "V", "ACN": "m", "ACT": "H", "AGN": "r", "AGT": "D",
       "ANT": "w", "CGN": "s", "CGT": "B", "CNT": "y", "GNT": "k",
       "*ACG": "v", "*ACN": "m", "*ACT": "h", "*AGN": "r", "*AGT": "d",
       "*ANT": "w", "*CGN": "s", "*CGT": "b", "*CNT": "y", "*GNT": "k",
       "ACGN": "v", "ACGT": "N", "ACNT": "h", "AGNT": "d", "CGNT": "b",
       "*ACGN": "v", "*ACGT": "N", "*ACNT": "h", "*AGNT": "d", "*CGNT": "b",
       "*ACGNT": "N"}


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
    f = open('{}/{}.genseries'.format(tmp_file_dir, col_name), 'a+')
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
    if opt.genfile_type == 'phylib':
        genseries_file = open(opt.file + '.gen_series.phy', 'w')
    elif opt.genfile_type == 'fasta':
        genseries_file = open(opt.file + '.gen_series.fasta', 'w')
    else:
        genseries_file = open(opt.file + '.gen_series.txt', 'w')

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

        # 同时创建缓存文件夹
        print(print_prefix + '缓存文件夹为:{}'.format(tmp_file_dir))
        if (not os.path.exists(tmp_file_dir)):
            os.mkdir(tmp_file_dir)

        # 考虑到内存太小，这里使用缓存文件作为中间处理
        for colname in sample_list:
            f = open('{}/{}.genseries'.format(tmp_file_dir, colname), 'w')
            f.close()

        # 处理进度相关
        progressbar = ProgressBar(maxval=opt.csv_size)
        progressbar.start()
        processed_size = 0

        # 启动进程池
        pool = mp.Pool(opt.parallel_jobs)

        # 按块读取csv文件，同时输出到对应基因序列中
        csv_reader = pd.read_csv(opt.file + '.csv', iterator=True, chunksize=opt.with_chunksize)
        for ck in csv_reader:
            # 估计当前的进度
            processed_size += ck.size * 2
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
                    f = open('{}/{}.genseries'.format(tmp_file_dir, sample_name), 'a+')
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
            sX_f = open('{}/{}.genseries'.format(tmp_file_dir, colname), 'r')

            if (opt.genfile_type == 'phylib'):
                genseries_file.write('seq{}\t'.format(i + 1))
            elif opt.genfile_type == 'fasta':
                genseries_file.write('>{}\n'.format(colname))
            genseries_file.write('{}\n'.format(sX_f.readline()))

        progressbar.finish()
        print(print_prefix + "文件合并完成")

        # 删除缓存文件夹，清理空间
        print(print_prefix + '清除缓存文件夹: del {}'.format(tmp_file_dir))
        shutil.rmtree(tmp_file_dir)

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

            # 针对不同格式的需求，写入tag信息
            if (opt.genfile_type == 'phylib'):
                genseries_file.write('seq{}\t'.format(i + 1))
            elif opt.genfile_type == 'fasta':
                genseries_file.write('>{}\n'.format(col_name))

            for GT in col.values:
                genseries_file.write(GT)
            genseries_file.write('\n')

        progressbar.finish()

    print(print_prefix + '基因序列生成完毕，序列总数:{}'.format(opt.sample_size))


ignore_ALTs = 0


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

    # 数值和碱基的对应关系
    nuc = {
        '0': s_data[3],
        '.': 'N',
    }

    for n, gtype in enumerate(s_data[4].split(',')):
        nuc[str(n + 1)] = gtype.replace('-', '*')

    # 从9列到最后，进行遍历处理
    for i in range(9, len(s_data)):
        # 获取基因型
        gentype = [s_data[i][0], s_data[i][2]]

        GT = []
        for i in gentype:
            GT.append(nuc[i])

        # 完成基因组合编码的转换。注意，这里使用set去除GT中的重复项
        out = amb[''.join(sorted(set(GT)))]

        if len(out) is 2:
            print("a")

        row += ',' + out

    return row


def covert_vcf2csv(opt):
    '''
    根据参数设置，解析数据。并生成csv文件
    :param opt: 使用的参数
    :return:
    '''
    src = open(opt.file, 'r')
    out = open(opt.file + '.csv', 'w')

    pline = ParallelLine(n_jobs=opt.parallel_jobs, chunk_size=opt.with_chunksize, show_process_status=True)
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
                                                 '--ignore_ALTs 1 --with_chunksize 10000 --parallel_jobs 4 ')
    parser.add_argument('--file', type=str, default='sample.vcf', help="需要处理的文件")
    parser.add_argument('--ignore_ALTs', type=int, default='0', help="变异类型超过该值则忽略该行。0表示不忽略")
    parser.add_argument('--covert_vcf2csv', type=str2bool, default=True,
                        help="是否执行covert_vcf2csv部分。如果已经转化出csv文件，设置true可以节省大量时间。可选参数：True,False")
    parser.add_argument('--gen_series', type=str2bool, default=True,
                        help="是否生成基因序列，建议设置基因型缺失时的替代符号，默认用空格代替。可选参数：True,False")
    parser.add_argument('--with_limited_mem', type=str2bool, default=True,
                        help="内存的使用策略。True尽量少的利用内存(速度相对慢一些)，False尽量多的利用内存(速度相对快一些)。可选参数：True,False")
    parser.add_argument('--with_chunksize', type=int, default=1000,
                        help="设置chunksize，控制分块读取时的行数。面对超大文件时可以有效降低内存使用。设置过大会延长加载时间")
    parser.add_argument('--parallel_jobs', type=int, default=4,
                        help="并行操作使用的进程数，设置1代表不并行。由于磁盘读取速度的限制,建议这里不要设置超过8。"
                             "如果采用了NVME的SSD，请检查线程的CPU利用率，CPU利用率低时，增加chunksize的大小，CPU利用率100%，可以适当增加parallel_jobs")
    parser.add_argument('--genfile_type', type=str, default='phylib',
                        help="生成的基因序列文件格式。可选有'txt'、'phylib'、'fasta'")

    opt = parser.parse_args()

    ignore_ALTs = opt.ignore_ALTs
    tmp_file_dir = '{}-tmp'.format(opt.file)

    print(opt)
    if (opt.covert_vcf2csv):
        covert_vcf2csv(opt)

    if (opt.gen_series):
        covert_csv2phylib(opt)

#!/usr/bin/env python3
"""
多物种基因组特征统计脚本（数值保留2位小数版）
用法: python genome_stats_multi.py --dir ./ --output all_results.tsv
"""

import sys
import os
import argparse
import numpy as np
from collections import defaultdict
import gzip
import glob
import re

def open_file(filename):
    """自动处理gzip压缩文件"""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def calculate_contig_stats(fasta_file):
    """计算组装序列的基本统计"""
    lengths = []
    total_gc = 0
    total_bases = 0
    
    with open_file(fasta_file) as f:
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    contig_seq = ''.join(seq).upper()
                    length = len(contig_seq)
                    lengths.append(length)
                    gc = contig_seq.count('G') + contig_seq.count('C')
                    total_gc += gc
                    total_bases += length
                seq = []
            else:
                seq.append(line)
        # 处理最后一个序列
        if seq:
            contig_seq = ''.join(seq).upper()
            length = len(contig_seq)
            lengths.append(length)
            gc = contig_seq.count('G') + contig_seq.count('C')
            total_gc += gc
            total_bases += length
    
    lengths.sort(reverse=True)
    total_len = sum(lengths)
    gc_content = (total_gc / total_bases * 100) if total_bases > 0 else 0
    
    # 计算N50/N90
    n50 = calculate_n(lengths, 50)
    n90 = calculate_n(lengths, 90)
    
    return {
        'Number of contigs': len(lengths),
        'Size (bp)': total_len,
        'GC content (%)': round(gc_content, 2),
        'N50 (bp)': n50,
        'N90 (bp)': n90,
        'Max length (bp)': max(lengths) if lengths else 0,
        'Min length (bp)': min(lengths) if lengths else 0
    }

def calculate_n(lengths, percent):
    """计算N50/N90"""
    total = sum(lengths)
    target = total * percent / 100
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= target:
            return l
    return 0

def parse_gff_genes(gff_file):
    """从GFF/GTF文件中提取基因结构信息"""
    transcripts = {}
    genes = {}
    
    with open_file(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feat_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attrs = parse_attributes(parts[8])
            
            feat_id = attrs.get('ID', '')
            parent = attrs.get('Parent', '')
            
            # 处理基因级别
            if feat_type == 'gene':
                gene_id = feat_id
                if gene_id:
                    genes[gene_id] = {
                        'gene_length': end - start + 1,
                        'transcripts': [],
                        'strand': parts[6],
                        'chr': parts[0]
                    }
            
            # 处理转录本/mRNA
            elif feat_type in ['mRNA', 'transcript']:
                trans_id = feat_id
                gene_id = parent if parent else attrs.get('gene_id', '')
                
                if not gene_id and '.' in trans_id:
                    gene_id = trans_id.split('.')[0]
                
                if trans_id:
                    transcripts[trans_id] = {
                        'exons': [],
                        'gene_id': gene_id,
                        'trans_length': end - start + 1
                    }
                    if gene_id and gene_id in genes:
                        genes[gene_id]['transcripts'].append(trans_id)
                    elif gene_id:
                        genes[gene_id] = {
                            'gene_length': 0,
                            'transcripts': [trans_id],
                            'strand': parts[6],
                            'chr': parts[0]
                        }
            
            # 处理外显子
            elif feat_type == 'exon':
                trans_id = parent
                if trans_id and trans_id in transcripts:
                    transcripts[trans_id]['exons'].append((start, end))
                elif trans_id and trans_id in genes:
                    if trans_id not in transcripts:
                        transcripts[trans_id] = {
                            'exons': [(start, end)],
                            'gene_id': trans_id,
                            'trans_length': end - start + 1
                        }
                        if trans_id not in genes:
                            genes[trans_id] = {
                                'gene_length': 0,
                                'transcripts': [trans_id],
                                'strand': parts[6],
                                'chr': parts[0]
                            }
            
            # 处理CDS
            elif feat_type == 'CDS':
                trans_id = parent
                if trans_id and trans_id in transcripts:
                    if (start, end) not in transcripts[trans_id]['exons']:
                        transcripts[trans_id]['exons'].append((start, end))
    
    # 如果没有找到任何转录本，尝试备用方法
    if not transcripts:
        transcripts, genes = parse_gff_alternative(gff_file)
    
    # 计算统计量
    gene_lengths = []
    exon_counts = []
    exon_lengths = []
    intron_lengths = []
    intron_counts = []
    
    for gene_id, gene_data in genes.items():
        trans_list = gene_data.get('transcripts', [])
        if not trans_list:
            continue
        
        # 选择最长的转录本
        best_trans_id = None
        best_trans_len = 0
        best_exons = []
        
        for trans_id in trans_list:
            if trans_id in transcripts:
                exons = transcripts[trans_id].get('exons', [])
                if exons:
                    trans_len = sum(e[1] - e[0] + 1 for e in exons)
                    if trans_len > best_trans_len:
                        best_trans_len = trans_len
                        best_trans_id = trans_id
                        best_exons = sorted(exons, key=lambda x: x[0])
        
        if not best_exons:
            continue
        
        # 基因长度
        gene_len = gene_data.get('gene_length', 0)
        if gene_len == 0 and best_exons:
            gene_len = best_exons[-1][1] - best_exons[0][0] + 1
        
        if gene_len > 0:
            gene_lengths.append(gene_len)
        
        # 外显子统计
        n_exons = len(best_exons)
        exon_counts.append(n_exons)
        for start, end in best_exons:
            exon_lengths.append(end - start + 1)
        
        # 内含子统计
        if n_exons > 1:
            n_introns = n_exons - 1
            intron_counts.append(n_introns)
            for i in range(n_introns):
                intron_start = best_exons[i][1] + 1
                intron_end = best_exons[i+1][0] - 1
                if intron_end >= intron_start:
                    intron_lengths.append(intron_end - intron_start + 1)
    
    return {
        'Average gene length (bp)': np.mean(gene_lengths) if gene_lengths else 0,
        'Average number of introns per gene': np.mean(intron_counts) if intron_counts else 0,
        'Average intron length (bp)': np.mean(intron_lengths) if intron_lengths else 0,
        'Average exons per gene': np.mean(exon_counts) if exon_counts else 0,
        'Average exon length (bp)': np.mean(exon_lengths) if exon_lengths else 0,
        'Number of genes': len(gene_lengths)
    }

def parse_gff_alternative(gff_file):
    """备用解析方法"""
    transcripts = {}
    genes = {}
    
    with open_file(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9 or parts[2] != 'exon':
                continue
            
            start = int(parts[3])
            end = int(parts[4])
            attrs = parse_attributes(parts[8])
            
            trans_id = attrs.get('transcript_id', attrs.get('Parent', attrs.get('ID', '')))
            gene_id = attrs.get('gene_id', trans_id.split('.')[0] if '.' in trans_id else trans_id)
            
            if trans_id:
                if trans_id not in transcripts:
                    transcripts[trans_id] = {'exons': [], 'gene_id': gene_id}
                transcripts[trans_id]['exons'].append((start, end))
                
                if gene_id not in genes:
                    genes[gene_id] = {'transcripts': [], 'gene_length': 0}
                if trans_id not in genes[gene_id]['transcripts']:
                    genes[gene_id]['transcripts'].append(trans_id)
    
    return transcripts, genes

def parse_attributes(attr_str):
    """解析GFF的attributes列"""
    attrs = {}
    for item in attr_str.split(';'):
        item = item.strip()
        if not item:
            continue
        if '=' in item:
            key, val = item.split('=', 1)
            val = val.strip('"')
            attrs[key] = val
        elif ' ' in item:
            parts = item.split(' ', 1)
            if len(parts) == 2:
                attrs[parts[0]] = parts[1].strip('"')
    return attrs

def process_single_species(species_name, genome_file, gff_file, verbose=True):
    """处理单个物种"""
    if verbose:
        print(f"\n处理物种: {species_name}")
        print(f"  基因组: {genome_file}")
        print(f"  注释: {gff_file}")
    
    # 检查文件是否存在
    if not os.path.exists(genome_file):
        print(f"  错误: 基因组文件不存在 - {genome_file}")
        return None
    if not os.path.exists(gff_file):
        print(f"  错误: 注释文件不存在 - {gff_file}")
        return None
    
    # 统计
    if verbose:
        print("  统计基因组组装指标...")
    contig_stats = calculate_contig_stats(genome_file)
    
    if verbose:
        print("  统计基因结构指标...")
    gene_stats = parse_gff_genes(gff_file)
    
    # 合并结果
    all_stats = {**contig_stats, **gene_stats}
    all_stats['Species'] = species_name
    
    # 重新排列列顺序，并格式化数值（保留2位小数）
    ordered_stats = {
        'Species': all_stats['Species'],
        'Number of contigs': all_stats['Number of contigs'],
        'Size (bp)': f"{all_stats['Size (bp)']:.2f}",
        'GC content (%)': f"{all_stats['GC content (%)']:.2f}",
        'N50 (bp)': f"{all_stats['N50 (bp)']:.2f}",
        'N90 (bp)': f"{all_stats['N90 (bp)']:.2f}",
        'Max length (bp)': f"{all_stats['Max length (bp)']:.2f}",
        'Min length (bp)': f"{all_stats['Min length (bp)']:.2f}",
        'Number of genes': all_stats['Number of genes'],
        'Average gene length (bp)': f"{all_stats['Average gene length (bp)']:.2f}",
        'Average exons per gene': f"{all_stats['Average exons per gene']:.2f}",
        'Average exon length (bp)': f"{all_stats['Average exon length (bp)']:.2f}",
        'Average number of introns per gene': f"{all_stats['Average number of introns per gene']:.2f}",
        'Average intron length (bp)': f"{all_stats['Average intron length (bp)']:.2f}"
    }
    
    if verbose:
        print(f"  完成! 发现 {all_stats['Number of contigs']} 个contigs, {all_stats['Number of genes']} 个基因")
    
    return ordered_stats

def find_paired_files(directory):
    """在目录中自动查找配对的基因组和注释文件（改进版）"""
    species_dict = {}
    
    # 获取所有文件
    all_files = os.listdir(directory)
    
    # 定义可能的基因组文件扩展名
    genome_extensions = ['.genome.fa', '.genome.fasta', '.fa', '.fasta', '.fna']
    gff_extensions = ['.gff', '.gff3', '.gtf']
    
    # 第一步：识别所有基因组文件
    genome_files = []
    for f in all_files:
        file_path = os.path.join(directory, f)
        if os.path.isfile(file_path):
            # 检查是否是基因组文件
            is_genome = False
            for ext in genome_extensions:
                if f.endswith(ext):
                    species_name = f.replace(ext, '')
                    genome_files.append((species_name, file_path))
                    is_genome = True
                    break
            # 也检查 .genome.fa.gz 这样的压缩文件
            if not is_genome and f.endswith('.genome.fa.gz'):
                species_name = f.replace('.genome.fa.gz', '')
                genome_files.append((species_name, file_path))
                is_genome = True
            if not is_genome and f.endswith('.genome.fasta.gz'):
                species_name = f.replace('.genome.fasta.gz', '')
                genome_files.append((species_name, file_path))
                is_genome = True
    
    # 第二步：识别所有GFF文件
    gff_files = []
    for f in all_files:
        file_path = os.path.join(directory, f)
        if os.path.isfile(file_path):
            is_gff = False
            for ext in gff_extensions:
                if f.endswith(ext):
                    species_name = f.replace(ext, '')
                    gff_files.append((species_name, file_path))
                    is_gff = True
                    break
            # 处理压缩文件
            if not is_gff and f.endswith('.gff.gz'):
                species_name = f.replace('.gff.gz', '')
                gff_files.append((species_name, file_path))
                is_gff = True
            if not is_gff and f.endswith('.gff3.gz'):
                species_name = f.replace('.gff3.gz', '')
                gff_files.append((species_name, file_path))
                is_gff = True
            if not is_gff and f.endswith('.gtf.gz'):
                species_name = f.replace('.gtf.gz', '')
                gff_files.append((species_name, file_path))
                is_gff = True
    
    # 第三步：匹配基因组和GFF文件
    paired_species = []
    
    print(f"找到 {len(genome_files)} 个基因组文件")
    print(f"找到 {len(gff_files)} 个注释文件")
    
    # 先尝试精确匹配物种名
    for genome_name, genome_path in genome_files:
        # 查找相同物种名的GFF文件
        matched_gff = None
        for gff_name, gff_path in gff_files:
            if genome_name == gff_name:
                matched_gff = (gff_name, gff_path)
                break
        
        # 如果没找到精确匹配，尝试部分匹配
        if not matched_gff:
            for gff_name, gff_path in gff_files:
                # 去除版本号等后缀后比较
                genome_base = re.sub(r'\.v\d+$', '', genome_name)
                gff_base = re.sub(r'\.v\d+$', '', gff_name)
                if genome_base == gff_base or genome_name in gff_name or gff_name in genome_name:
                    matched_gff = (gff_name, gff_path)
                    break
        
        if matched_gff:
            paired_species.append((genome_name, genome_path, matched_gff[1]))
            print(f"  匹配成功: {genome_name} -> {os.path.basename(genome_path)} + {os.path.basename(matched_gff[1])}")
        else:
            print(f"  警告: 未找到 {genome_name} 对应的GFF文件")
    
    return paired_species

def read_batch_config(config_file):
    """读取批量配置文件"""
    species_list = []
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                species_name = parts[0]
                genome_file = parts[1]
                gff_file = parts[2]
                species_list.append((species_name, genome_file, gff_file))
            elif len(parts) == 2:
                # 假设格式: genome.fasta gff.gff，使用文件名作为物种名
                species_name = os.path.splitext(os.path.basename(parts[0]))[0]
                genome_file = parts[0]
                gff_file = parts[1]
                species_list.append((species_name, genome_file, gff_file))
    return species_list

def main():
    parser = argparse.ArgumentParser(description='多物种基因组特征统计')
    parser.add_argument('--single', nargs=3, metavar=('SPECIES', 'GENOME', 'GFF'),
                        help='单个物种模式: 物种名 基因组.fasta 注释.gff')
    parser.add_argument('--batch', metavar='CONFIG_FILE',
                        help='批量模式: 配置文件（每行: 物种名\t基因组.fasta\t注释.gff）')
    parser.add_argument('--dir', metavar='DIRECTORY',
                        help='目录模式: 自动查找目录中的配对文件')
    parser.add_argument('--output', '-o', default='genome_stats.tsv',
                        help='输出文件（默认: genome_stats.tsv）')
    parser.add_argument('--verbose', '-v', action='store_true', default=True,
                        help='显示详细信息')
    
    args = parser.parse_args()
    
    species_data = []
    
    # 单物种模式
    if args.single:
        species_name, genome_file, gff_file = args.single
        species_data = [(species_name, genome_file, gff_file)]
    
    # 批量模式
    elif args.batch:
        if not os.path.exists(args.batch):
            print(f"错误: 配置文件不存在 - {args.batch}")
            sys.exit(1)
        species_data = read_batch_config(args.batch)
        print(f"从配置文件加载了 {len(species_data)} 个物种")
    
    # 目录模式
    elif args.dir:
        if not os.path.exists(args.dir):
            print(f"错误: 目录不存在 - {args.dir}")
            sys.exit(1)
        print(f"正在扫描目录: {args.dir}")
        species_data = find_paired_files(args.dir)
        print(f"在目录中找到了 {len(species_data)} 个物种")
    
    else:
        print("错误: 请指定 --single, --batch 或 --dir 模式")
        parser.print_help()
        sys.exit(1)
    
    if not species_data:
        print("错误: 没有找到任何物种数据")
        sys.exit(1)
    
    # 处理所有物种
    results = []
    for species_name, genome_file, gff_file in species_data:
        result = process_single_species(species_name, genome_file, gff_file, args.verbose)
        if result:
            results.append(result)
    
    if not results:
        print("错误: 没有成功处理任何物种")
        sys.exit(1)
    
    # 写入输出文件
    if results:
        # 获取所有列名
        all_columns = list(results[0].keys())
        
        with open(args.output, 'w') as f:
            # 写入表头
            f.write('\t'.join(all_columns) + '\n')
            
            # 写入数据
            for result in results:
                row = [str(result.get(col, 'NA')) for col in all_columns]
                f.write('\t'.join(row) + '\n')
        
        print(f"\n所有结果已保存到: {args.output}")
        
        # 打印汇总（数值保留2位小数）
        print("\n===== 统计汇总 =====")
        print(f"{'Species':25} {'Contigs':>8} {'Size(Mb)':>10} {'GC%':>8} {'N50(kb)':>10} {'Genes':>8}")
        print("-" * 75)
        for result in results:
            # 从格式化的字符串中提取数值用于显示
            size_mb = float(result['Size (bp)']) / 1_000_000
            n50_kb = float(result['N50 (bp)']) / 1000
            gc_pct = float(result['GC content (%)'])
            print(f"{result['Species']:25} {result['Number of contigs']:>8} {size_mb:>10.2f} "
                  f"{gc_pct:>8.2f} {n50_kb:>10.2f} {result['Number of genes']:>8}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
合并多个基因计数文件的脚本
用法: python merge_gene_counts.py file1.txt file2.txt [file3.txt ...] -o output.txt
"""

import argparse
import pandas as pd
from pathlib import Path

def read_gene_file(filepath):
    """读取基因计数文件，返回字典"""
    df = pd.read_csv(filepath, sep='\t')
    # 假设第一列是基因名，第二列是计数值
    gene_col = df.columns[0]
    count_col = df.columns[1]
    
    # 从文件名提取样本名（去掉.txt后缀）
    sample_name = Path(filepath).stem
    
    return df.set_index(gene_col)[count_col].to_dict(), sample_name

def merge_gene_files(file_list, output_file, fill_value=0):
    """合并多个基因文件"""
    all_data = {}
    sample_names = []
    
    # 读取所有文件
    for filepath in file_list:
        gene_dict, sample_name = read_gene_file(filepath)
        all_data[sample_name] = gene_dict
        sample_names.append(sample_name)
    
    # 获取所有唯一的基因名
    all_genes = set()
    for gene_dict in all_data.values():
        all_genes.update(gene_dict.keys())
    
    # 创建合并的数据框
    merged_df = pd.DataFrame(index=sorted(all_genes))
    
    # 填充每个样本的数据
    for sample_name in sample_names:
        merged_df[sample_name] = merged_df.index.map(
            lambda x: all_data[sample_name].get(x, fill_value)
        )
    
    # 添加行总和（可选）
    merged_df['Total'] = merged_df.sum(axis=1)
    
    # 保存结果
    merged_df.to_csv(output_file, sep='\t')
    print(f"合并完成！共处理 {len(file_list)} 个文件，{len(all_genes)} 个基因")
    print(f"结果已保存到: {output_file}")
    
    return merged_df

def main():
    parser = argparse.ArgumentParser(description='合并多个基因计数文件')
    parser.add_argument('files', nargs='+', help='输入的基因计数文件')
    parser.add_argument('-o', '--output', default='merged_genes.txt', 
                       help='输出文件名 (默认: merged_genes.txt)')
    parser.add_argument('-f', '--fill', type=int, default=0,
                       help='缺失基因的填充值 (默认: 0)')
    
    args = parser.parse_args()
    
    merge_gene_files(args.files, args.output, args.fill)

if __name__ == '__main__':
    main()

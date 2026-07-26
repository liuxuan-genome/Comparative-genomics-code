#!/usr/bin/env python3
"""
批量处理所有.cazy.txt文件
"""

from pathlib import Path

def batch_add_headers():
    """批量处理所有CAZy文件"""
    cazy_files = list(Path('.').glob('*.cazy.txt'))
    
    if not cazy_files:
        print("未找到 .cazy.txt 文件")
        return
    
    print(f"找到 {len(cazy_files)} 个文件")
    
    for filepath in cazy_files:
        # 提取物种名
        species_name = filepath.stem.split('.')[0]
        
        # 读取内容
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # 检查是否已有表头
        if lines and 'Gene' in lines[0]:
            print(f"跳过 {filepath.name} (已有表头)")
            continue
        
        # 添加表头
        with open(filepath, 'w') as f:
            f.write(f"Gene\t{species_name}\n")
            f.writelines(lines)
        
        print(f"✓ 处理完成: {filepath.name} -> 表头: Gene\t{species_name}")

if __name__ == '__main__':
    batch_add_headers()

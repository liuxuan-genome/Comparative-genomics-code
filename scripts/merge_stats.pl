#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use File::Spec;

# 用法：perl merge_stats.pl

my @result_dirs = glob("*.result");  # 查找所有 .result 目录
my %data;
my @categories = qw(GT CE GH AA CBM PL);
my @species_list;

foreach my $dir (@result_dirs) {
    # 从目录名提取物种名（去掉 .result 后缀）
    my $species = $dir;
    $species =~ s/\.result$//;
    push @species_list, $species;
    
    # 构建统计文件路径
    my $stats_file = File::Spec->catfile($dir, "cazy_statistics.txt");
    
    next unless -f $stats_file;  # 如果文件不存在则跳过
    
    open my $fh, "<", $stats_file or die "Cannot open $stats_file: $!";
    
    # 读取文件内容
    while (my $line = <$fh>) {
        chomp $line;
        
        # 找到数据行（包含冒号的）
        if ($line =~ /^([A-Z]+)\s*:\s*(\d+)/) {
            my $cat = $1;
            my $count = $2;
            if (grep { $_ eq $cat } @categories) {
                $data{$species}{$cat} = $count;
            }
        }
        # 找到总计行
        if ($line =~ /^总计\s*:\s*(\d+)/) {
            $data{$species}{'Total'} = $1;
        }
    }
    
    close $fh;
    
    # 补全缺失的分类（计数为0）
    foreach my $cat (@categories) {
        $data{$species}{$cat} //= 0;
    }
    $data{$species}{'Total'} //= 0;
}

# 输出合并后的表格
print "Species";
foreach my $cat (@categories) {
    print "\t$cat";
}
print "\tTotal\n";

foreach my $species (@species_list) {
    print $species;
    foreach my $cat (@categories) {
        print "\t" . ($data{$species}{$cat} // 0);
    }
    print "\t" . ($data{$species}{'Total'} // 0);
    print "\n";
}

# 输出所有物种合计
print "\n" . "=" x 80 . "\n";
print "所有物种合计:\n";
print "=" x 80 . "\n";
my %total_all;
foreach my $species (@species_list) {
    foreach my $cat (@categories) {
        $total_all{$cat} += $data{$species}{$cat};
    }
    $total_all{'Total'} += $data{$species}{'Total'};
}

print "Category";
foreach my $cat (@categories) {
    print "\t$cat";
}
print "\tTotal\n";
print "All";
foreach my $cat (@categories) {
    print "\t" . ($total_all{$cat} // 0);
}
print "\t" . ($total_all{'Total'} // 0);
print "\n";

import igraph as ig
import leidenalg
import pandas as pd
import time
# 记录程序开始时间
start_time = time.time()
# 读取 k-plex 社区
kplex_communities = []
with open('kplex_communities.txt', 'r') as f:
    for line in f:
        community = list(map(int, line.strip().split()))
        kplex_communities.append(community)

# print("K-plex communities:",len(kplex_communities))
# print("K-plex communities:", kplex_communities)


#读取零散节点
scattered_nodes = []
for i in range(34):
    scattered_nodes.append(i)

#print(scattered_nodes)
# with open('scattered_node.txt', 'r') as f:
#     for line in f:
#         scattered_nodes.append(int(line.strip()))

# print("Scattered nodes:", scattered_nodes)

# 定义文件路径
file_path = './Dataset/zarachy.txt'  # 替换为实际文件路径
#file_path = './Dataset/football.txt'  # 替换为实际文件路径
#file_path = './Dataset/email-Eu-core.txt'  # 替换为实际文件路径
#file_path = './Dataset/dblp-small.txt'  # 替换为实际文件路径
#file_path = './Dataset/dblp-large.txt'  # 替换为实际文件路径

# 初始化 edge_list
edge_list = []

# 读取文件并构建边列表
with open(file_path, 'r') as file:
    for line in file:
        nodes = line.strip().split()  # 去除行末尾的换行符并分割
        if len(nodes) == 2:  # 确保每行有两个节点
            u, v = map(int, nodes)  # 将节点转换为整数
            edge_list.append((u, v))  # 将边存储到 edge_list 中

# # 输出读取的边列表
# print("Edge list:", edge_list)


# 构建 igraph 图对象
g = ig.Graph(edges=edge_list, directed=False)

# initial_membership = [-1] * len(kplex_communities)

# 初始社区划分：k-plex 社区
initial_membership = [-1] * g.vcount()  # 初始化所有节点的社区为 -1
community_id = 0

# 为每个 k-plex 社区分配社区 ID
for community in kplex_communities:
    for node in community:
        initial_membership[node] = community_id
    community_id += 1

# 将零散节点分配给自己的单独社区
for node in scattered_nodes:
    initial_membership[node] = community_id
    community_id += 1



# 使用 Leiden 算法优化社区
partition = leidenalg.find_partition(g, leidenalg.ModularityVertexPartition, initial_membership=initial_membership)

# 记录程序结束时间
end_time = time.time()

# 输出程序的运行时间
print(f"运行时间: {end_time - start_time} 秒")


# # 输出优化后的社区结果
# for comm_id, community in enumerate(partition):
#     print(f"Community {comm_id}: {community}")



# 初始的社区划分 (k-plex) 已存储在 initial_membership 中
# Leiden 算法优化后的社区划分 (partition) 需要转换成一个向量

# 创建 Leiden 算法后的社区标签
final_membership = [-1] * g.vcount()

for comm_id, community in enumerate(partition):
    for node in community:
        final_membership[node] = comm_id


# 定义处理函数，读取文件并提取第二个数值
def extract_second_column(file_path):
    second_column_list = []
    
    # 打开并读取文件
    with open(file_path, 'r') as file:
        for line in file:
            # 只提取每行的第二个数值
            _, second_value = map(int, line.strip().split())
            second_column_list.append(second_value)
    
    return second_column_list

# 示例文件路径，假设你的文件为 "community_data.txt"
file_path = './data_labels/zarachy-labels.txt'
#file_path = './data_labels/football-labels.txt'
#file_path = './data_labels/email-Eu-core-department-labels.txt'
#file_path = './data_labels/dblp-small-labels.txt'
#file_path = './data_labels/dblp-large-labels.txt'

# 处理文件并获取第二列的数据
labels_true  = extract_second_column(file_path)


from sklearn.metrics import normalized_mutual_info_score
# 计算 NMI
nmi_score = normalized_mutual_info_score(labels_true, final_membership)

# 输出 NMI
print(f"NMI score: {nmi_score}")

from sklearn.metrics import f1_score
from collections import defaultdict
import numpy as np

# 计算每个预测标签对应的真实标签
label_mapping = defaultdict(lambda: defaultdict(int))

# 统计每个预测标签对应的真实标签的计数
for true, pred in zip(labels_true, final_membership):
    label_mapping[pred][true] += 1

# # 打印标签映射
# for pred_label, true_labels in label_mapping.items():
#     print(f"Predicted Label {pred_label}: {dict(true_labels)}")

# 根据计数决定最可能的真实标签
best_mapping = {}
for pred_label, true_labels in label_mapping.items():
    best_true_label = max(true_labels, key=true_labels.get)
    best_mapping[pred_label] = best_true_label

# 使用映射来替换 final_membership 中的标签
mapped_final_membership = [best_mapping[pred] for pred in final_membership]

# 计算 F1-score
f1 = f1_score(labels_true, mapped_final_membership, average='weighted')

# 输出 F1-score
print(f"F1-score after mapping: {f1}")



def jaccard_similarity(setA, setB):
    """计算两个集合的 Jaccard 相似度"""
    intersection = len(set(setA) & set(setB))  # 交集大小
    union = len(set(setA) | set(setB))          # 并集大小
    return intersection / union if union != 0 else 0  # 返回 Jaccard 相似度

# 将真实标签和预测标签按标签分组
true_communities = defaultdict(set)
predicted_communities = defaultdict(set)

for node, label in enumerate(labels_true):
    true_communities[label].add(node)

for node, label in enumerate(mapped_final_membership):         #  ???是否对应标签
    predicted_communities[label].add(node)

# 计算最终的 Jaccard 相似度
total_intersection = 0
total_union = 0

# 计算每个社区的 Jaccard 相似度并累加
for pred_label, pred_nodes in predicted_communities.items():
    if pred_label in true_communities:  # 确保预测标签在真实标签中
        true_nodes = true_communities[pred_label]
        total_intersection += len(pred_nodes & true_nodes)  # 交集
        total_union += len(pred_nodes | true_nodes)          # 并集

# 计算最终的 Jaccard 相似度
final_jaccard_score = total_intersection / total_union if total_union != 0 else 0

# 输出最终的 Jaccard 相似度
print(f"Final Jaccard Similarity: {final_jaccard_score}")

# 打开并读取文件
with open('Dataset/zarachy.txt', 'r') as file:
    lines = file.readlines()

# 初始化空的边列表
edges = []

# 遍历每行，找到所有的 edge
for i in range(len(lines)):
    if 'source' in lines[i]:
        source = int(lines[i].split()[1])  # 提取 source 的值
        target = int(lines[i+1].split()[1])  # 提取 target 的值
        edges.append((source, target))  # 将边加入边列表

# 输出边列表，以所需的格式
with open('Dataset/zarachy.txt', 'w') as out_file:
    for edge in edges:
        out_file.write(f"{edge[0]-1} {edge[1]-1}\n")

print("finish!")

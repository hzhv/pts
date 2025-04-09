import matplotlib.pyplot as plt

def load_csr_matrix(file_path):
    """
    从文件中读取CSR矩阵数据，文件格式要求按行存储：
      1. n
      2. row_ptr.size()
      3. row_ptr (所有元素用空格隔开)
      4. col_id.size()
      5. col_id (所有元素用空格隔开)
      6. val.size()
      7. val (所有元素用空格隔开)
    """
    with open(file_path, 'r') as f:
        n = int(f.readline().strip())
        
        rp_size = int(f.readline().strip())
        # 假设所有row_ptr数据位于一行中（数值间以空格隔开）
        row_ptr_line = f.readline().strip()
        row_ptr = list(map(int, row_ptr_line.split()))
        if len(row_ptr) != rp_size:
            raise ValueError("row_ptr的元素数目与声明的不符！")
        
        col_size = int(f.readline().strip())
        col_id_line = f.readline().strip()
        col_id = list(map(int, col_id_line.split()))
        if len(col_id) != col_size:
            raise ValueError("col_id的元素数目与声明的不符！")
            
        val_size = int(f.readline().strip())
        val_line = f.readline().strip()
        val = list(map(float, val_line.split()))
        if len(val) != val_size:
            raise ValueError("val的元素数目与声明的不符！")
            
    return {"n": n, "row_ptr": row_ptr, "col_id": col_id, "val": val}

def level_scheduling(csr):
    """
    以CSR格式的下三角矩阵，执行层调度算法：
      - 只考虑 j < i 且 A(i,j) ≠ 0 的情况
      - 返回一个二维列表，每个子列表表示一层中包含的节点/行号
    """
    n = csr["n"]
    row_ptr = csr["row_ptr"]
    col_id = csr["col_id"]
    
    # 计算每行的入度；依赖关系：对于每个非零 A(i,j) 且 j < i，认为 i依赖于 j
    in_degree = [0] * n
    dependents = [[] for _ in range(n)]
    
    for i in range(n):
        # 遍历行 i 的所有非零项：下标范围 row_ptr[i] 到 row_ptr[i+1]-1
        for idx in range(row_ptr[i], row_ptr[i+1]):
            j = col_id[idx]
            # 仅处理下三角（j < i）的依赖关系
            if j < i:
                in_degree[i] += 1
                dependents[j].append(i)
    
    levels = []
    current_level = [i for i in range(n) if in_degree[i] == 0]
    
    # 分层调度：不断删除当前层节点，对其所有依赖节点减入度
    while current_level:
        levels.append(current_level)
        next_level = []
        for node in current_level:
            for dep in dependents[node]:
                in_degree[dep] -= 1
                if in_degree[dep] == 0:
                    next_level.append(dep)
        current_level = next_level
    
    # # 检查是否存在环（如果有节点入度不为 0，则依赖图中存在环）
    # if any(degree > 0 for degree in in_degree):
    #     raise RuntimeError("检测到依赖图中存在环！")
        
    return levels

def plot_histogram(levels, filename="level_histogram"):
    """
    根据分层结果绘制直方图：
      - 横坐标：层的索引
      - 纵坐标：每一层中包含的行数
    返回每层行数的列表。
    """
    level_counts = [len(level) for level in levels]
    
    top_n = 200
    plt.figure(figsize=(10, 6))
    plt.bar(range(top_n), level_counts[:top_n])
    plt.xlabel("Level Index")
    plt.ylabel("Number of Rows")
    plt.title("Top 1000 Levels (Most Active)")
    # plt.tight_layout()

    plt.savefig(filename + ".png")
    
    return level_counts

def main():
    file_path = "L_FEM_3D_thermal2_csr.txt"
    file_name = file_path.split("/")[-1].replace(".txt", "")
    print("Loading Large Sparse A in CSR format...")
    csr = load_csr_matrix(file_path)
    print("Loaded. Start level scheduling...")
    
    levels = level_scheduling(csr)
    total_levels = len(levels)
    print(f"Level scheduling done. Found {total_levels} levels in total.")

    level_counts = plot_histogram(levels, file_name)
    
    average_rows = sum(level_counts) / total_levels if total_levels > 0 else 0
    print(f"Each level contains {average_rows:.2f} rows on average.")    

if __name__ == "__main__":
    main()

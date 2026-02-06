import numpy as np
import matplotlib.pyplot as plt
import os

# 配置
FILENAME = "build/samples.bin"  # 确保路径正确
N_DIM = 512
SAMPLE_DIM = 2 * N_DIM  # s1 和 s2 拼接
DTYPE = np.int16

def verify_covariance():
    if not os.path.exists(FILENAME):
        print(f"Error: File {FILENAME} not found. Run ./dump_samples first.")
        return

    print("Loading data...")
    # 读取二进制文件
    raw_data = np.fromfile(FILENAME, dtype=DTYPE)
    
    # 重塑数组: (Num_Samples, 1024)
    num_samples = raw_data.size // SAMPLE_DIM
    print(f"Loaded {num_samples} samples. Dimension: {SAMPLE_DIM}")
    
    X = raw_data.reshape((num_samples, SAMPLE_DIM))

    # Step 1: 估计样本协方差矩阵
    # Sigma_hat = Cov(X)
    print("Computing Covariance Matrix...")
    Sigma_hat = np.cov(X, rowvar=False)

    # Step 2: 估计方差 sigma^2
    # 理论上应该是对角矩阵，对角线元素即为方差
    sigma2_hat = np.mean(np.diag(Sigma_hat))
    sigma_hat = np.sqrt(sigma2_hat)
    
    print("-" * 30)
    print(f"Estimated Variance (sigma^2): {sigma2_hat:.4f}")
    print(f"Estimated Std Dev (sigma) : {sigma_hat:.4f}")
    print("-" * 30)

    # Step 3: 构造理想矩阵 Sigma_ideal = sigma^2 * I
    Sigma_ideal = sigma2_hat * np.eye(SAMPLE_DIM)

    # Step 4: 计算误差 (Frobenius Norm)
    # 归一化误差，除以元素个数以便观察
    diff = Sigma_hat - Sigma_ideal
    error_norm = np.linalg.norm(diff, ord='fro')
    avg_error = error_norm / (SAMPLE_DIM) 

    print(f"Frobenius Norm of Difference: {error_norm:.4f}")
    print(f"Average Error per element   : {avg_error:.6f}")

    # 判定标准
    # 理想情况下，非对角元素应该接近 0，对角元素接近 sigma^2
    off_diag = Sigma_hat.copy()
    np.fill_diagonal(off_diag, 0)
    max_off_diag = np.max(np.abs(off_diag))
    
    print(f"Max Off-Diagonal Correlation: {max_off_diag:.4f}")

    # Step 5: 可视化
    plt.figure(figsize=(10, 8))
    plt.title(f"Covariance Matrix Heatmap (N={num_samples})\n$\\sigma^2 \\approx {sigma2_hat:.1f}$")
    
    # 使用 imshow 显示矩阵
    # 理想情况应该是一条明亮的对角线，背景全黑
    im = plt.imshow(Sigma_hat, cmap='viridis', interpolation='nearest')
    plt.colorbar(im, label='Covariance')
    plt.xlabel("Component Index (0-511: s1, 512-1023: s2)")
    plt.ylabel("Component Index")
    
    plt.tight_layout()
    plt.show()

    # 额外：画个直方图看看分布形状
    plt.figure(figsize=(10, 6))
    plt.hist(X.flatten(), bins=100, density=True, alpha=0.6, color='g', label='Sample Data')
    
    # 拟合高斯曲线对比
    x_range = np.linspace(X.min(), X.max(), 1000)
    pdf = (1/(sigma_hat * np.sqrt(2 * np.pi))) * np.exp(-0.5 * (x_range / sigma_hat)**2)
    plt.plot(x_range, pdf, 'r--', linewidth=2, label=r'Ideal Gaussian ($\sigma={sigma_hat:.2f}$)')
    
    plt.title("Distribution of Coefficients (s1, s2)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()

if __name__ == "__main__":
    verify_covariance()
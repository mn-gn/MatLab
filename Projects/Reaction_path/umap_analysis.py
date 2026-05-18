import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import umap
from mpl_toolkits.mplot3d import Axes3D

# 1. 파일 경로 설정 및 데이터 로드
# 이전 MATLAB 스크립트에서 저장한 .mat 파일의 절대 경로를 입력합니다.
mat_filepath = r"C:\Users\ming0\Desktop\MATLAB\Reaction_path\UMAP_Input_Data.mat"

if not os.path.exists(mat_filepath):
    print(f"파일을 찾을 수 없습니다: {mat_filepath}")
    exit()

print("데이터를 로드하는 중입니다...")
data = sio.loadmat(mat_filepath)

# MATLAB에서 저장한 변수들 추출 (.flatten()을 사용하여 1차원 배열로 변환)
time_vector = data['time_vector'].flatten()
umap_input_RSV = data['umap_input_RSV']
umap_input_SRSV = data['umap_input_SRSV']
num_components = data['num_components'].flatten()[0]

# 시간 벡터를 색상으로 매핑하기 위해 로그 스케일 적용
# (음수 딜레이나 0 근처의 값이 있을 수 있으므로 절대값과 작은 상수를 더함)
color_var = np.arange(len(time_vector))

# 2. UMAP 분석 설정 및 실행
# 반응 중간체의 군집화 정도를 조절하는 파라미터
n_neighbors_val = 50
min_dist_val = 0.1

print("UMAP 임베딩을 계산 중입니다. 잠시만 기다려주세요...")

# RSV (Right Singular Vectors) UMAP 계산
reducer_rsv_2d = umap.UMAP(n_components=2, n_neighbors=n_neighbors_val, min_dist=min_dist_val, random_state=42)
umap_2d_rsv = reducer_rsv_2d.fit_transform(umap_input_RSV)

reducer_rsv_3d = umap.UMAP(n_components=3, n_neighbors=n_neighbors_val, min_dist=min_dist_val, random_state=42)
umap_3d_rsv = reducer_rsv_3d.fit_transform(umap_input_RSV)

# SRSV (Scaled Right Singular Vectors) UMAP 계산
reducer_srsv_2d = umap.UMAP(n_components=2, n_neighbors=n_neighbors_val, min_dist=min_dist_val, random_state=42)
umap_2d_srsv = reducer_srsv_2d.fit_transform(umap_input_SRSV)

reducer_srsv_3d = umap.UMAP(n_components=3, n_neighbors=n_neighbors_val, min_dist=min_dist_val, random_state=42)
umap_3d_srsv = reducer_srsv_3d.fit_transform(umap_input_SRSV)

# 3. 2D 및 3D 비교 시각화 (2x2 서브플롯)
fig = plt.figure(figsize=(16, 12))
fig.canvas.manager.set_window_title('UMAP Trajectories: RSV vs SRSV')

# [1] RSV - 2D
ax1 = fig.add_subplot(2, 2, 1)
scatter1 = ax1.scatter(umap_2d_rsv[:, 0], umap_2d_rsv[:, 1], c=color_var, cmap='jet', s=40, alpha=0.8)
fig.colorbar(scatter1, ax=ax1, label='log10(Time)')
ax1.set_title('RSV - 2D UMAP')
ax1.set_xlabel('UMAP 1')
ax1.set_ylabel('UMAP 2')

# [2] RSV - 3D
ax2 = fig.add_subplot(2, 2, 2, projection='3d')
scatter2 = ax2.scatter(umap_3d_rsv[:, 0], umap_3d_rsv[:, 1], umap_3d_rsv[:, 2], c=color_var, cmap='jet', s=40, alpha=0.8)
fig.colorbar(scatter2, ax=ax2, label='log10(Time)', pad=0.1)
ax2.set_title('RSV - 3D UMAP')
ax2.set_xlabel('UMAP 1')
ax2.set_ylabel('UMAP 2')
ax2.set_zlabel('UMAP 3')

# [3] SRSV - 2D
ax3 = fig.add_subplot(2, 2, 3)
scatter3 = ax3.scatter(umap_2d_srsv[:, 0], umap_2d_srsv[:, 1], c=color_var, cmap='jet', s=40, alpha=0.8)
fig.colorbar(scatter3, ax=ax3, label='log10(Time)')
ax3.set_title('SRSV - 2D UMAP')
ax3.set_xlabel('UMAP 1')
ax3.set_ylabel('UMAP 2')

# [4] SRSV - 3D
ax4 = fig.add_subplot(2, 2, 4, projection='3d')
scatter4 = ax4.scatter(umap_3d_srsv[:, 0], umap_3d_srsv[:, 1], umap_3d_srsv[:, 2], c=color_var, cmap='jet', s=40, alpha=0.8)
fig.colorbar(scatter4, ax=ax4, label='log10(Time)', pad=0.1)
ax4.set_title('SRSV - 3D UMAP')
ax4.set_xlabel('UMAP 1')
ax4.set_ylabel('UMAP 2')
ax4.set_zlabel('UMAP 3')

plt.suptitle(f'UMAP Trajectories Over Time (Top {num_components} Components)', fontsize=16)
plt.tight_layout()
plt.show()
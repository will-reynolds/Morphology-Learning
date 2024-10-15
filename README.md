# Cell-growth-on-scaffold-by-morphology-learning 
#### Our project is composed of two parts: morphology learning and structural optimization.

![Overview (1)](https://github.com/user-attachments/assets/baf17543-c8aa-42c6-a355-c0a3032afb93)
[Weiming Wang](https://weiming-wang.github.io/), [Yanhao Hou](https://scholar.google.com/citations?user=hFebGnoAAAAJ&hl=zh-CN), [Renbo Su](https://www.researchgate.net/profile/Renbo_Su), [Weiguang Wang](https://scholar.google.co.uk/citations?user=z60XpXQAAAAJ&hl=en), [Charlie C.L.Wang](https://mewangcl.github.io/). Simultaneously optimized mechanical stiffness and cell growth on scaffold by morphology learning.


## Abstract
The morphology of microstructures in scaffolds is a critical factor influencing both mechanical properties and biological performance, such as cell growth, in bone tissue engineering. While metamaterials with optimized mechanical properties, such as hyperelasticity, thermal efficiency, or energy absorption, have been designed through complex microstructural geometry with the help of numerical simulation of multi-physics, optimizing biological performance alongside these physical traits remains a significant challenge. This is primarily because biological performance must be evaluated through laboratory experiments. To address this issue, we propose a novel approach that simultaneously optimizes mechanical stiffness and cell growth effectiveness through a data-driven morphology learning method. Our approach effectively extracts and learns shape patterns from a specifically selected dataset of microstructures recognized for their superior cell growth performance. We then integrate these patterns into a computational framework based on topology optimization, enabling the generation of new microstructures that maintain high levels of cell growth while enhancing mechanical stiffness. Using this morphology learning-based method, we identified a scaffold design for bone tissue engineering that demonstrates a significant improvement in both mechanical stiffness (by 26.35\%) and cell proliferation.

## Requirements
#### Numpy 1.18.1
#### Pytorch 1.5.0
#### scipy 1.4.1
#### cvxopt 1.2.0

## Usage

# Step 1: Download dataset and sample points:

Before that, you should download the full dataset or the filtered dataset from the following links:

[Full dataset](https://drive.google.com/file/d/13bQZYW0j-nxYzCS1y7LEDq5f7ejf5fNt/view?usp=sharing)

[Filtered dataset](https://drive.google.com/file/d/1bQKuHTiDrEkyLNP7c-Ne0qq0PsQvRqMh/view?usp=sharing)

Then, a set of points as well as their SDF values should be sampled for each structure in the dataset. 
You can run the following file to sample points:

Run ./code/Morphology Learning/data/sampling

# Step 2: Morphology learning:

You can run the following file to train the neural network:

Run ./code/Morphology Learning/Morphology Learning.py

# Step 3: Structural optimization:

After training, you can run the following file to optimize structural performance:

Run ./code/Topology Optimization/structural_opt.py

## Citations:


## Contact information:

Weiming Wang (weiming.wang@manchester.ac.uk)

Yanhao Hou (yanhao.hou@manchester.ac.uk)

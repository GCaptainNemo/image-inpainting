# 基于结构稀疏的图像修复算法
## 简介
图像修复(image inpainting)是一种将受损图像恢复至图像原始状态的技术。非数据驱动的图像修复算法可以理解为利用图像完好区域对图像损坏区
域进行推断的过程。具体来说就是首先对待修复图像像素进行一定的分布假设，然后构造统计量对损坏的像素进行点估计。以均值滤波来说，可以认为图像以某像素为中心的一个邻域内，服从以该像素值为均值的一个未知分布，利用旁边完整的像素构造的统计量    T(X) = (x1 + x2 + x3 + ... + XN)/N     对受损像素进行估计。
理论上所有非数据驱动的图像修复算法都可以统一到如上的框架中，但有些算法并不能显式地写出假设分布的形式，这里介绍利用结构稀疏(structure sparsity)作为先验进行图像修复。


## 特点
    基于结构稀疏的图像修复算法对受损直线边界的修复有较好的表现，算法以patch为单位进行修复可以较完美地复原纹理(texture)和结构(structure)。对于非直
线边界，该方法是无能为力的；此外，如果直线边界过长且过程中存在光照、纹理的渐变，该方法并不能捕捉到该渐变，修复后可能会导致一定视觉上不自然。这时只能使
用数据驱动的修复算法，比如深度学习、字典学习等等。

## 效果
![results](./results/horse_results.png)


## 参考文献
[1] Jin D , Bai X . Patch-Sparsity-Based Image Inpainting Through a Facet Deduced Directional Derivative[J]. IEEE Transactions on Circuits
and Systems for Video Technology, 2019, 29(5):1310-1324.

[2] Z. Xu and J. Sun, “Image inpainting by patch propagation using patch sparsity,” IEEE Trans. Image Process., vol. 19, no. 5, pp. 1153–1165,
May 2010.

[3] A. Criminisi, P. Perez, and K. Toyama, “Region filling and object removal by exemplar-based image inpainting,” IEEE Trans. Image Process., 
vol. 13, no. 9, pp. 1200–1212, Sep. 2004. 

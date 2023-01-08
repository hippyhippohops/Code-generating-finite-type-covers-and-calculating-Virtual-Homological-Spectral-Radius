# Hyperbolic-Volume-and-Virtual-Homological-Spectral-Radius-Conjecture

Let $\Sigma_{g,n}$ be an orientable surface of genus $g$ and $n$ punctures, such that $2-2g-n < 0$. Let $f \in $ Mod($\Sigma_{g,n}$) denote an element in the mapping class group. The study of the bounds of $D(f)$ and $D_{h}(f)$ is a dynamic field in geometric topology. \\

Given any automorphisms $f$ of a connected orientable surface $\Sigma_{g,n}$, it is evident that any virtual homological spectral radius for $f$ is greater than or equal to $1$. It is shown by C.T. McMullen that any virtual homological spectral radius for a pseudo-Anosov surface automorphism $f$ is strictly lesser than the dilation if the invariant foliations for $f$ have prong singularities of odd order \cite{McMullen}. Since the topology of the mapping torus depends only on the mapping class $f$, we denotes its topological type by $N_{f}$. A celebrated theorem by Thurston \cite{Thurston} asserts that $N_{f}$ admits a hyperbolic structure iff $f$ is pseudo-Anosov. Kojiwa-Macshane proved that:
\begin{align*}
    ln(D(f)) \geq \frac{Vol(M_f)}{3 \pi |\chi(X)|},
\end{align*}
where $M_f$, $\chi(X)$ and $D(f)$ are the mapping torus, Euler characteristic of the surface and the dilation of $f$ respectively \cite{Kojima-Mcshane}. \\


With references to these, Dr Thang Le conjectured relations between the spectral radius of $\widetilde{f}$ acting on $H_{1}(\widetilde{\Sigma_{g,n}},\mathbb{Z})$, and the volume of the mapping torus with respect to $f$, 
\begin{align*}
    M_f := \frac{\mathbb{I} \text{ x } X}{(1,x) \sim (0,f(x))},
\end{align*}
where $\widetilde{\Sigma_{g,n}}$ varies over the finite covers of $\Sigma_{g,n}$ to which $f$ lifts. \\

For my senior thesis project, I developed this code to construct all possible finite type covers of the figure $8$, given an integer $n$. The program will then obtain the spectral radius from each finite type cover from the collection of all covers constructed and test the conjectures. All these theory will be developed in this thesis, together with a brief introduction to Mapping Class Groups and Pseudo-Anosov Maps. Then, I will explain how this theory is translated to work on the code. Results of the $4$-sheeted cover will be explained thoroughly step-by-step and the results for $5$ and $6$ sheeted covers will also be included. Finally, I will explain the checks we have put in place to ensure the correctness of the algorithm before I conclude this thesis with suggestions for improvement. \\

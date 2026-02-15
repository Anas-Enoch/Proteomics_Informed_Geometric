\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{figure1.pdf}
\caption{
\textbf{Overview of the metabolite-centric geometric framework for metabolic networks.}
(A) A schematic bipartite representation of a metabolic network, with metabolites and reactions connected by stoichiometric relationships. This representation highlights that metabolites are coupled indirectly through shared reactions rather than through direct adjacency.
(B) Condition-specific reaction coupling strengths inferred from gene-level proteomics via gene–protein–reaction (GPR) rules. Differences between conditions modulate reaction weights without altering network topology.
(C) Induced geometric embedding of metabolites derived from the low-frequency spectrum of the metabolite Laplacian, illustrating how mechanistic connectivity defines distances between metabolites.
(D) Partial metabolomic observability, in which only a subset of metabolites is measured. Observed metabolites (Ω) define a restricted operator whose ability to preserve the full geometry motivates the identifiability criterion used throughout the study.
}
\label{fig:overview}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{figure2.pdf}
\caption{
\textbf{Proteomics-informed operator construction and geometric interpretation.}
(A) Construction of the condition-specific metabolite Laplacian
$\Delta_M^{(c)} = W_M^{1/2} S W_R^{(c)} S^\top W_M^{1/2}$,
combining stoichiometric structure with reaction coupling weights derived from
proteomics and metabolite reliability weighting.
(B) Block structure of the associated Dirac operator $D^{(c)}$, which couples
metabolite and reaction spaces through the weighted stoichiometric operator
$d^{(c)}$.
(C) Squaring the Dirac operator yields block-diagonal Laplacians on metabolites
and reactions, making explicit that metabolite geometry emerges from
stoichiometry-mediated reaction coupling rather than from direct metabolite
adjacency.
(D) Intuition linking coupling strength to geometric distance: metabolites
connected through strongly weighted reactions are geometrically close, while
weakly coupled metabolites are separated by larger diffusion distances.
}
\label{fig:operator}
\end{figure}

\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{figure3.pdf}
\caption{
\textbf{Geometry preservation under partial observability.}
(A) Error in global spectral signatures as a function of metabolite panel size for random selection, topology-based baselines, and the proposed geometry-aware framework.
(B) Distribution of pairwise diffusion distance errors under partial observability, showing improved preservation of local geometry by the proposed method.
(C) Comparison of geometry preservation with and without proteomics-informed reaction weighting, demonstrating that condition-specific proteomics improves robustness to missing metabolite measurements.
(D) Example pathway-level geometry illustrating how the restricted operator preserves relative metabolite organization when informative panels are selected.
}
\label{fig:identifiability}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{figure4.pdf}
\caption{
\textbf{Geometry-aware metabolite panel design.}
(A) Greedy selection procedure that iteratively adds metabolites to minimize geometric distortion under partial observability.
(B) Comparison of panel efficiency across selection strategies, showing that geometry-aware selection achieves lower distortion with fewer metabolites.
(C) Trade-off between robust panels shared across conditions and condition-adaptive panels optimized for specific proteomic states.
(D) Pathway membership of selected metabolites, illustrating that geometry-preserving panels preferentially sample mechanistically central pathways rather than high-degree metabolites.
}
\label{fig:panel}
\end{figure}


\begin{figure}[H]
\centering
\includegraphics[width=\linewidth]{Fig5_permutation_robustness.pdf}
\caption{
\textbf{Biological alignment of proteomics-driven reaction weights 
stabilizes metabolic geometry under partial observation.}

(A) Geometry preservation versus panel size under different selection strategies.
(B) Sensitivity envelope across continuous OR aggregators.
(C) Robustness across nonlinear saturation mappings.
(D) Disease classification performance using diffusion embeddings.
(E) Permutation null distribution demonstrating collapse of geometric stability 
when proteomic weights are randomly reassigned.
}
\label{fig:Biological alignment}
\end{figure}
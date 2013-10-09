Description
=================================================

An en-route monitoring agency (EMA) monitors the safe operation of Performance-based Navigation (PBN) aircraft in the international airspace. Collision risk models are often used for the periodic safety assessment by EMA. Lateral deviations of aircraft from the intended route are modelled as a probability distribution. Gaussian distributions, Laplace distributions and the mixture of them are often used as the distribution model. The models are estimated from the lateral deviation data observed by air surveillance systems such as radars and mandate reports for large deviations. Aircraft flying in the oceanic airspace are allowed to conduct 1 nautical mile or 2 nautical miles intentional lateral offsets without informing air traffic controls. When we develop the above distribution model, we should take into account these facts. The difficulty in the data analysis is that there are no means to know whether intentional offset is applied or unintentionally deviated when an aircraft deviates 1 nautical mile right from the assigned route. We developed the parameter estimation algorithm for the lateral deviation model simultaneously estimating the applied intentional lateral offset. Please check the following paper for the details.

M. Fujita, Estimation of Navigation Performance and Offset by the EM algorithm and the variational Bayesian method, Advances and Applications in Statistics, vol. 35, no. 1, pp.1--27, 2013.

This program provides a Java-based code implementing the parameter estimation algorithm described above.


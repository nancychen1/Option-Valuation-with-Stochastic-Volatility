# Option Valuation with Stochastic Volatility
 Values call options via Monte Carlo simulation using a GARCH stochastic volatility model

Using an initial stock price of 100, and continuously compounding risk free rate of 5% yearly, with time till expiration of 0.5 years, the GARCH model is used to generate the stock price paths. The long run volatility is 30% annually and stochastic volatility starts off at 35% annually. 
Using the model with three cases:
1.  alpha = beta = 0.0, gamma = 1.0
2.  alpha = .01, beta = 0.0 and gamma = 0.99
3.  alpha = 0.01, beta = .10, and gamma = 0.89

We value the call options with strikes of 60 to 180 in steps of 10 and then using Black-Scholes, we compute a call option's implied volatility. Then, we plot the results of the strike vs the implied volatility. 

We see in the first two cases, we have constant or near constant implied volatility, the volatility is deterministic, and in the third case we see the shape of the graph mimics the real world volatility skew/smile. This Volatility Smile shows that as an option becomes more Out of Money (OTM) or In the Money (ITM), the underlying asset's implied volatility increases. The stochastic nature of real-world volatility is a contributing factor in the smile.   

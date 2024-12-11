In this write-up for all the surprises I encountered trying to be clever with the making change problem, I will describe the techniques I learned along the way, including two surprises: generating functions, and convolutions using the fast Fourier transform.

# The problem
Given some monetary amount X, can you count how many unique ways there are to combine bills and coins that add up to X?

## Dynamic programming
Change-making is a quintessential [dynamic programming problem](https://en.wikipedia.org/wiki/Change-making_problem) right alongside [matrix chain multiplication](https://en.wikipedia.org/wiki/Matrix_chain_multiplication) and any problem where the solution to the entire problem can be related to solutions of smaller sub-problems (called [optimal substructure](https://en.wikipedia.org/wiki/Optimal_substructure)). Unlike [divide-and-conquer](https://en.wikipedia.org/wiki/Divide-and-conquer_algorithm), the sub-problems of a dynamic programming solution overlap. I appreciate the examples Wikipedia uses to describe this dichotomy: quicksort and merge sort are divide-and-conquer while Djikstra's shortest path and recursively calculating a Fibonacci sequence are dynamic programming.

### Memo-y
Given the overlap, dynamic programming uses [memoization](https://en.wikipedia.org/wiki/Memoization) so every substep is only ever computed once and then cached, and combinatorial expansion avoided. Using the example of calculating the *nth* Fibonacci number, memoization reduces the [complexity](https://en.wikipedia.org/wiki/Computational_complexity_theory) of `Fib(n)` from $O(\phi^n)$ to just $O(n)$.

Of course this is wonderful, what a vast improvement right?

## Hot take
Well, sure... but in practice there's overhead to memoization. There's a cost, and it's not *just* memory: a shred of clarity is lost. Granted, this is a very minor contrivance, though I can't help but shed a tear for the elegant recursion like the classic 1970's littering PSA.

Another "problem" with that vast improvement, (for the specific case of the *nth* Fibonacci number) is that there's a simple, non-DP, non-recursive, generalized expression for calculating the value, known as [Binet's formula](https://en.wikipedia.org/wiki/Fibonacci_sequence#Relation_to_the_golden_ratio) (this is where that $\phi$ came from in the analysis above). Just two irrational constants are needed to achieve $O(1)$ complexity.

There is an argument to be made that the recursive approach avoids loss of precision, however, the unfortunate truth is that this too is bested by an even more fool-proof method: the humble lookup table. For all the clever approaches to calculating a Fibonacci number, there's none faster, more accurate, or more straightforward. The $O(\phi^n)$ growth means that any method used would quickly reach the limit of even 64-bit unsigned integers; a lookup table needs only 93 entries (assuming the zeroth entry starts at 0) with the final being $[92]=7,540,113,804,746,346,429$.

While it just so happens that the Fibonacci lookup table corellates 1 to 1 with the dynamic approach, this is possibly the only case where that is true.

### Analytic solutions
Any mathematician or computer scientist can appreciate the joy of taking a $O(...)$ problem and deriving a $O(1)$ solution. I've heard this called an "analytic solution" before which seems appropriate, however it's not the exact same as [the mathematical definition of the same-named class of expressions](https://en.wikipedia.org/wiki/Closed-form_expression#Comparison_of_different_classes_of_expressions) since this definition necessarily involves a complexity analysis (mathematically speaking, the recursive solution is a single finite sum of ones for any given input, making it an algebraic expression). Nevertheless, I use the name here for the sake of brevity.

While not every problem has an analytic solution, when one exists it is often the best, fastest, canonically preferred means of solving such a problem. Having experience with the dynamic programming solution to the change-making problem, I'd always guessed there could be a more analytic solution to the problem, some sort of formula where you plug in the monetary amount and it just returns an answer never fussing with recursion or memoization. This intuition comes from a plot of the dynamic programming function.

<img src="https://github.com/user-attachments/assets/ccf20e12-aaae-4f55-aafe-2afe2fe092b2" height="900" />

You'll quickly notice (for USD) that the number of ways to make change only increases on $0.05 increments, and further still, it looks like there *could* be some sort of pattern there too for the dime and quarter and so on... The numbers grow with each denomination in surprising ways, and yes **there is an exploitable pattern**.

# The analysis
Though I'd guessed at it initially, I didn't attempt to solve it analytically before stumbling across [this YouTube video from Mathologer](https://youtu.be/VLbePGBOVeg). In it professor Burkard Polster opens with an unexpected insight that leads to a clever specialized solution.

The unexpected insight? The YouTube algorithm is scary good, but so are [generating functions](https://en.wikipedia.org/wiki/Generating_function). I had first seen generating functions in [a different video by 3Blue1Brown](https://youtu.be/bOXCLR3Wric), where he uses slightly different techniques to solve a different counting problem, and that's worth a watch too.

I will explain the gist, but the Mathologer video explains everything step-by-step in exquisite detail should you wish to see it.

## Generating functions
Suppose you have three polynomials, for example

$$P(x) = 1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^{10}$$
$$N(x) = 1 + x^5 + x^{10}$$
$$D(x) = 1 + x^{10}$$

Now suppose that I suggest you multiply them all together, you'll get something like

$$USC(x) = P(x)N(x)D(x) = $$
$$1 + x + x^2 + x^3 + x^4 + 2x^5 + 2x^6 + 2x^7 + 2x^8 + 2x^9 + 4x^{10} +$$
$$3x^{11} + 3x^{12} + 3x^{13} + 3x^{14} + 4x^{15} + 3x^{16} + 3x^{17} + 3x^{18} + 3x^{19} + 4x^{20} +$$
$$2x^{21} + 2x^{22} + 2x^{23} + 2x^{24} + 2x^{25} + x^{26} + x^{27} + x^{28} + x^{29} + x^{30}$$

If you consider the coefficient for $x^{10}$, namely $4$, you might notice that there are *4* distinct ways to give change for $0.10 in $USC$: 10 pennies, 5 pennies and 1 nickel, 2 nickels, and finally 1 dime. This is no coincidence, as it corresponds exactly to how polynomial products work.

In this sense, the coefficient for a given power of $x$ in any polynomial product is an answer to the question "how many unique ways can the terms of the starting polynomials combine to achieve this power of $x$ (scaled by initial coefficients)?" This actually genuinely surprised me to learn, but **polynomial multiplication is a form of [convolution](https://en.wikipedia.org/wiki/Convolution)**. This will come up again later.

This is guaranteed to work for all powers of $x$ up to the degree of those starting polynomials. "Neat," I imagine you thinking quietly to yourself if you're unfamiliar with generating functions. "Is there a trivially easy way to get answers for numbers larger than $0.10 USD?" you might ask and the answer would be yes and then no.

The astute reader would notice the equations above were constructed to model individual denominations, and it is trivially easy to extend them. Just add more and more powers of $x$ at the appropriate intervals for their respective denomination up to the power representing the value you're after, multiply and combine like terms, find the coefficient for the appropriate power of $x$, and Bob's your uncle. Solved.

### "Wait, that's terrible as an algorithm"
"Straight polynomial multiplication is $O(d_1d_2)$ on the degrees of the polynomials. That solution requires multiplying a series of D (count of different denominations) polynomials each potentially the size of N (the input value), giving a complexity of $O(N^D)$. That's not even analytic as you defined it."

I like the cut of your jib, hypothetical reader. It's true there's work yet to turn an intractable solution into something usable, daresay even tractable.

First, we can formalize those polynomials given above. Each denomination in our hypothetical currency was expressed as a sum of powers of $x$ up to the amount to give change for. To give change for any possible amount, we can create infinite versions of the sums (called [series](https://en.wikipedia.org/wiki/Series_(mathematics))) of the form

$$\sum_{i=0}^{\infty} x^{id} = 1 + x^d + x^{2d} + x^{3d} + ...$$

where $d$ is the denomination in the smallest units (cents, fen, paisa, etc). These infinite series are called [generating functions](https://en.wikipedia.org/wiki/Generating_function) when synthesized purely for the behavior of polynomial operations. Now, with some clever algebra and a whiff of calculus we can use the formula for an [infinite geometric sum](https://en.wikipedia.org/wiki/Geometric_series#Convergence_of_the_series_and_its_proof) to express each of the generating functions as a simple fraction. If you've never seen the formula for an infinite geometric sum, it is

$$\sum_{i=0}^{\infty} x^i = 1 + x + x^2 + x^3 + ... = \frac{1}{1-x}$$

which only works for values $-1 < x < 1$. This is of little concern however since we will never need to *evaluate* these polynomials. Applying this formula to the product of our generating functions (for the six denominations of US coinage) comes to

$$\frac{1}{1-x}\frac{1}{1-x^5}\frac{1}{1-x^{10}}\frac{1}{1-x^{25}}\frac{1}{1-x^{50}}\frac{1}{1-x^{100}}$$

While this looks more manageable, this form is only really useful as-is if we were trying to evaluate this polynomial for inputs of our synthetic variable $x$. In order for it to be of any use to find coefficients, we need to convert it back to a polynomial form. The first step would be to combine like terms.

### The unit coin
If you noticed that all the powers of $x$ except the first are some multiple of $5$, well spotted. This is an optimization that made sense for me to implement in code for reasons that will become very obvious in a bit. First, let's rewrite that first polynomial in terms of $x^5$.

$$\Big(\frac{1+x+x^2+x^3+x^4}{1+x+x^2+x^3+x^4}\Big)\Big(\frac{1}{1-x}\Big)=(1+x+x^2+x^3+x^4)\frac{1}{1-x^5}$$

What this component of our polynomial product really represents is whether or not our particular currency has a 'unit coin'. That $(1+x+x^2+x^3+x^4)$ portion of this product captures the only means of achieving an exponent that is not a multiple of $5$. This corresponds exactly with the idea that if our currency doesn't have a unit coin, (the Canadian dollar for instance), that there would be no way to combine individual components of the generating functions to achieve a power of $x$ of $1$.

### The rest of the denominations
For now, we will tuck that 'unit coin polynomial' away for now to focus on the *remainder* of the products

$$\Big(\frac{1}{1-x^5}\Big)^2\frac{1}{1-x^{10}}\frac{1}{1-x^{25}}\frac{1}{1-x^{50}}\frac{1}{1-x^{100}}$$

For each one, the very same trick can be used to convert each term to $\frac{1}{1-x^{100}}$ such as for the quarter

$$\frac{1+x^{25}+x^{50}+x^{75}}{1+x^{25}+x^{50}+x^{75}}\frac{1}{1-x^{25}}=(1+x^{25}+x^{50}+x^{75})\frac{1}{1-x^{100}}$$

Doing this for all of them will result in

$$(1+x^5+x^{10}+...+x^{95})^2(1+x^{10}+x^{20}+...+x^{90})(1+x^{25}+x^{50}+x^{75})(1+x^{50})\frac{1}{(1-x^{100})^6}$$

Next we can take all those extracted polynomial numerators and expand them into their own 'remainder polynomial' like so

$$1+2x^5+4x^{10}+6x^{15}+9x^{20}+13x^{25}+...+980x^{195}+985x^{200}+985x^{205}+980x^{210}+...+13x^{380}+9x^{385}+6x^{390}+4x^{395}+2x^{400}+x^{405}$$

You *could* multiply in the 'unit coin polynomial' now, but the effect would be to "stretch" the 'remainder polynomial' out, duplicating its coefficients by a factor of the number of terms of the 'unit coin polynomial'. **It is in this way that extracting that 5 term unit coin polynomial lets us contract a 406 term polynomial into an 82 term polynomial.** This is a necessary optimization as we will see later.

### Generalized geometric sum
There is also a more generalized form of the geometric sum above which uses [binomial coefficient notation](https://en.wikipedia.org/wiki/Binomial_coefficient), a version raised to a power $n$ like this

$$\sum_{i=0}^{\infty}{{n-1+i}\choose{n-1}}x^i={{n-1}\choose{n-1}}+{{n}\choose{n-1}}x+{{n+1}\choose{n-1}}x^2+...=\frac{1}{(1-x)^n}$$

again only for $-1<x<1$. This formula looks like a nice fit for the polynomial product above if we swap $x$ for $x^{100}$ and $n=6$. For our initial case of US coinage, we finally, *finally* end with the product of three things, the finite 5 term 'unit coin polynomial', the finite 81 term 'remainder polynomial', and a new infinite polynomial sum that altogether looks like this

$$(1+x+x^2+x^3+x^4)\times$$
$$(1+2x^5+4x^{10}+6x^{15}+9x^{20}+13x^{25}+18x^{30}+...+x^{405})\times$$
$$[{5\choose5}+{6\choose5}x^{100}+{7\choose5}x^{200}+{8\choose5}x^{300}+{9\choose5}x^{400}+...]$$

## How to use it
With the above, we need to only consider the values of coefficients that add up to a power of $x$ that matches our desired input monetary amount. For example, the video gives $4.00 USD, or 400 cents, so this is equivalent to finding all possible combinations of terms that multiply to $x^{400}$. There are many relevant terms, and here I will isolate them

$$(1+...)\times$$
$$(1+...+287x^{100}+...+985x^{200}+...325x^{300}+...2x^{400}+...)\times$$
$$[{5\choose5}+{6\choose5}x^{100}+{7\choose5}x^{200}+{8\choose5}x^{300}+{9\choose5}x^{400}+...]$$

For all possible input amounts, there is always one and only one relevant term in the first polynomial, so it can be completely ignored for now. Notice how the answer can be given by multiplying and combining the remaining visible terms

$$1{9\choose5}x^{400}+287x^{100}{8\choose5}x^{300}+985x^{200}{7\choose5}x^{200}+325x^{300}{6\choose5}x^{100}+{5\choose5}2x^{400}$$

the variables all combine to $x^{400}$ except it isn't relevant anymore so we can remove it for clarity

$${9\choose5}+287{8\choose5}+985{7\choose5}+325{6\choose5}+{5\choose5}=38,835$$

and therefore there are 38,835 ways to give change for $4.00 USD. That's it! Having a finite polynomial up front allows us to limit the search for coefficients in the infinite polynomial. For *any* non-negative integer amount, to find the number of ways to give change, you only need to find all the ways to combine these terms of the infinite polynomial **over the range of the finite polynomial**, and it extends in a very manageable way.

The Mathologer video shows that the formula for whole dollar amounts greater than 3 is

$$\\#(k\$)={k+5\choose5}+287{k+4\choose5}+985{k+3\choose5}+325{k+2\choose5}+{k+1\choose5}$$

In code if we combined the two finite polynomials it would look like:

```cpp
uint64_t analyticWaysToGiveChange(uint64_t amount) const {
    uint64_t count = 0ul;
    for (uint64_t exponent = amount % 100ul; exponent < coefficients.size() && exponent <= amount; exponent += 100ul)
        count += coefficients[exponent] * nCk(((amount - exponent) / 100ul) + 5ul, 5ul);

    return count;
}
```

In my code above, `nCk(n, k)` is ${n \choose k}$. This works for any amount, but eagle-eyed reeders would notice that there are a lot of magic numbers in there and that I misspelled 'readers'. There's a few more loose ends to tie up before reaching a final algorithm.

## Loose ends

### The skeleton in the closet
Remember earlier when a hypothetical reader pointed out that "straight polynomial multiplication is $O(d_1d_2)$"? These are sparse polynomials, meaning they contain low information, but multiplying is still not trivial. Some shortcuts can be taken, but ultimately, it's just a much longer form of FOIL.

This means that the setup for building the finite portion of the polynomial is not so innocent as it might seem.

### Naming constants
In the snippet above, I'll give a few names to these numbers. First, the `5` is just 1 less than the number of different denominations. There is a subtlety to that `100`, however, it's not merely the largest denomination. In fact it's a constant that's crucial to the construction of the initial finite portion of the polynomial that was simply emergent from the initial setup.

Imagine a version without the dollar or half dollar coins. In this setup, our generating functions will ultimately reach the form of

$$\frac{1}{1-x}\frac{1}{1-x^5}\frac{1}{1-x^{10}}\frac{1}{1-x^{25}}$$

and when we try to unify all these together the way we did before, we'll find that in order to combine them in the same manner as before, we have to convert each one to a denominator $1-x^{50}$.

$$\frac{1+x+...+x^{49}}{1-x^{50}}\frac{1+x^5+...+x^{45}}{1-x^{50}}\frac{1+x^{10}+...+x^{40}}{1-x^{50}}\frac{1+x^{25}}{1-x^{50}}$$

That exponent is the [least common multiple](https://en.wikipedia.org/wiki/Least_common_multiple) of all the different denominations. The bigger the LCM and the more denominations you have, the bigger the finite polynomial is. In fact, it can be said that the size of the finite polynomial is directly dependant on the number of denominations and how many unique factors they all share. The worst case is when all the denominations of a hypothetical currency are coprime to each other.

### The GCD
In the case of the Canadian and United States dollars, without their pennies, there would be no way to give exact change for amounts between 5 cent totals. However this is not true for the Euro, where with its mighty 2 cent and 5 cent coins there's only two denominations that cannot be achieved without its penny, €0.01 and €0.03.

This idea can be captured by finding the [greatest common divisor](https://en.wikipedia.org/wiki/Greatest_common_divisor) of all the denominations greater than $1$. If the GCD is still $1$ then there may be only a few gaps, however the GCD being greater than one leads to two major consequences. First, that's what gave the USD its 5 cent steps, where the counts can only differ between the $0.05 boundaries. The second follows from the first, and it is that this pattern can be exploited to reduce the size of the finite polynomial, making the multiplications faster.

I find this last fact somewhat amusing considering that this is a speedup employed earlier to make the math easier to do by hand. In similar situations usually it is a complication to employ such tricks as computers can often be faster than such optimizations through sheer brute force.

This is unfortunately not one of those times. To give a taste of some numbers, when constructing the finite polynomial for the USD, our polynomial can be contracted down to a manageable ~16,000 coefficients. The Euro, however, which cannot be contracted since its GCD is 1 cent, has more than 1.4 million coefficients in its finite polynomial.

## Why would you think to do this in order to solve this problem?
The Mathologer video credits Graham, Knuth, and Patashnik's [*Concrete Mathematics*](https://en.wikipedia.org/wiki/Concrete_Mathematics) for inspiring the use of generating functions and indeed chapter 7.1 is all about this *exact* application for generating functions. While all of those authors are accomplished mathematicians, the preface of the book explains the provenance for its material reaching even further back in time.

Even the aforementioned 3Blue1Brown video admits difficulty finding a way to motivate how someone might discover generating functions on their own. Personally, I think the spark of inspiration lies in the doing: in order to even have a chance to make connections one must first know of something to connect. Though I had seen that 3B1B video before, undertaking this solution truly galvanized them for me.

## The code
With all those changes above, we arrive at more or less the final version I use in my solution
```cpp
uint64_t analyticWaysToGiveChange(uint64_t amount) const {
    if (!m_hasUnitValue && amount % m_gcd != 0ul)
        return 0ul;

    const uint64_t k = m_denominations.size() - 1ul;
    const uint64_t limit = (coefficients.size()) * m_gcd;
    uint64_t count = 0ul;
    for (uint64_t exponent = amount % m_lcm; exponent < limit && exponent <= amount; exponent += m_lcm)
        count += m_coefficients[exponent / m_gcd] * nCk(((amount - exponent) / m_lcm) + k, k);

    return count;
}
```

In this version, at startup we split the unit value (if it exists) into the condensed form in the polynomial and save just a single boolean to handle it with a quick if statement up front. Next we pick the "choose" portion of our binomial coefficient, `k` which is the same for all terms. Now we need to know when to stop searching the infinite portion of the polynomial, this needs to happen when we run out of finite polynomial coefficients, but since we're potentially using a condensed polynomial, we just multiply it by the GCD, the scaling factor. This factor shows up again to ensure we index the right coefficient. To wrap it all up, the LCM is the step for the exponents in the infinite polynomial, so this is the step in the loop.

## Initial results
Very, very promising... so long as you do not include the "startup" cost of multiplying the polynomials together. I was impressed by just how fast it was; so much that aggregate timing is required to measure the analytic solution, even for very large inputs that take the dynamic approach roughly four or five orders of magnitude longer.

My initial approach to multiplying the polynomials was laughably slow, however, and was quickly thrown out. I mentioned this to a few friends of mine, and got two very different answers. "Compute the polynomial at compile-time," said one which gave me a good chuckle. This is 100% a very reasonable solution to do for a real-world problem, but in this case, I wanted a means of doing this solution *faster* than a dynamic implementation, and precomputation would sort of cheat this goal. However, there is something to be said that this is a very natural extension that would be very hard to do with the dynamic programming approach.

Another friend suggested I use the fast Fourier transform and I was instantly sold.

# The FFT
The Fourier Transform is so interrelated with so many areas of mathematics that I do not think it's possible for the breadth of its versatility to sink in all at once. Such comprehension is a gradual process.

Though I've used the FFT before, this is the first time I've written one myself. This is what I would have liked to know before using the FFT.

## *Quick* disambiguation of topics from Fourier transforms for those unfamiliar

### Fourier series and analysis
A [Fourier series](https://en.wikipedia.org/wiki/Fourier_series) is the summation of signals, (sinusoids, frequencies, etc.) used in [Fourier analysis](https://en.wikipedia.org/wiki/Fourier_analysis) which is what most people likely think about when they think of Joseph Fourier. These infinite sums of sinusoidal signals can be truncated to approximate arbitrary functions, evaluating only finitely many. If you are familiar with [animations involving vectors aligned tip-to-tail spinning to trace out complex patterns](https://youtu.be/-qgreAUpPwM), these are accomplished using Fourier series approximations.

The series and its function approximations can be created using the Fourier transform.

### The (continuous) Fourier transform
The [Fourier transform](https://en.wikipedia.org/wiki/Fourier_transform) is a formula used to express a continuous 'input' function $f(x)$ as a continuous 'output' function $\widehat{f}(\xi)$ which describes the *frequencies* of the input.

The formula for the transform looks like

$$\widehat{f}(\xi)=\int_{-\infty}^\infty f(x)e^{-i2\pi\xi x}dx$$

The transform works by using [Euler's formula](https://en.wikipedia.org/wiki/Euler's_formula) to "wind" the input around the origin of the complex plane and finding the "center of mass" so to speak (the integral). That "center of mass" for a given $\xi$ (greek xi, the frequency of the winding) is $\widehat{f}(\xi)$. By changing $\xi$, the transform "winds" $f(x)$ tighter or more loosely and the "center of mass" moves, $\widehat{f}(\xi)$ changes. When $\widehat{f}(\xi)$ is very near the origin, $f(x)$ has low correlation with this particular "winding", $\xi$, and when the center is far from the origin, there is high correlation.

For an animated explanation of the above, I recommend [this 3Blue1Brown video](https://youtu.be/spUNpyF58BY) for FT novices like I was.

While the domain of $f(x)$ can be either $\mathbb{R}$ or $\mathbb{C}$, $\widehat{f}(\xi)$ is a complex output; the magnitude encodes the strength of the correlation, and the [argument](https://en.wikipedia.org/wiki/Argument_(complex_analysis)) of $\widehat{f}(\xi)$ encodes the phase. If your data are all in phase, however, it is common to discard the imaginary component.

### The discrete Fourier transform
The [Discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform) (conventionally denoted calligraphic F, $\mathcal{F}$) does the same thing as its continuous counterpart using discrete inputs (usually called samples) to extract discrete frequencies. The DFT is effectively the same formula as the FT but swapping an integral for a summation, and operating on discrete (usually *finite*) inputs instead of continuous functions. The number of discrete frequencies it can extract is equivalent to the number of input samples.

While the continuous Fourier transform is more commonly thought of as a **formula**, the DFT can be considered as both a formula and an **algorithm**. A DFT computes the sum of $N$ inputs "evaluated" at $N$ distinct points to yield $N$ outputs. A naive approach to "evaluate" $N$ inputs at $N$ points would be $O(N^2)$.

### Inverse DFT
This is a function that does what it says on the tin, in particular it satisfies the following: $\mathcal{F}^{-1}(\mathcal{F}(f))=f$ for input function $f$. When *computing* $\mathcal{F}$ and $\mathcal{F}^{-1}$, (perhaps surprisingly) there is very little distinction between the two (just one exponent). For this reason they are sometimes treated as the same operation, particularly in the context of computing. This remains true for the FFT.

### The FFT
The [fast Fourier transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) is a family of algorithms to calculate the DFT. While the FFT is still named after Joseph Fourier, the algorithms were first invented some time in the 19th or 20th century well after his death.

The simplest and most common version is [the Cooley-Tukey FFT](https://en.wikipedia.org/wiki/Cooley–Tukey_FFT_algorithm) which uses a divide-and-conquer approach to exploit symmetry by evaluating the inputs on the [roots of unity](https://en.wikipedia.org/wiki/Root_of_unity). What's clever is that **using the roots of unity as its evaluation points is what induces that symmetry**. The intuition for how this works can be seen in [this video from Reducible](https://youtu.be/h7apO7q16V0), but the gist is that it takes the $O(N^2)$ run time down to $O(N\log(N))$.

In practice, the FFT is a complete drop-in replacement for the DFT in many applications; it is rare to find a use-case where the general DFT is necessary over the FFT.

## The FFT and polynomials
I was always so used to the Fourier transform being the "signal processing" function that when I first learned that it could be used to multiply polynomials, I was a little surprised. Now for me, it is practically the definition of what the FFT is. To understand how, first we need to know two quirks about polynomial multiplication.

### Linearity
For two polynomials $f(x)$ and $g(x)$ and any input value $i$, multiplying the outputs $f(i)$ and $g(i)$ together is equivalent to evaluating on the combined multiplied polynomial, e.g. $f(x)\times g(x)=(f\times g)(x)$.

### Unique interpolation
Given $N+1$ input-output pairs in the form of $(x_n,f(x_n))$, there is unique degree-$N$ polynomial that defines $f(x)$. For example given the pairs $(-1, 3), (1, -3), (3, -1)$, there is one and only one degree-$2$ polynomial that passes through those points, $f(x)=x^2-3x-1$.

Deriving a polynomial from such a set of points is called [polynomial interpolation](https://en.wikipedia.org/wiki/Polynomial_interpolation) and there is a variety of means to this end.

### Putting it together
In the disambiguation above, I described the DFT as "computing" or "evaluating" a set of inputs on "points". This is *exactly* analogous to evaluating a polynomial, because that *is* what the DFT is doing.

If "evaluation" is something that takes a set of inputs and returns the outputs, then the "inverse" of this, something that takes those outputs and returns the associated inputs, this is by definition "interpolation". It then follows that the FFT "evaluates" and the IFFT "interpolates".

## Everything is polynomials
Multiplying two polynomials together perfectly models a [discrete convolution](https://en.wikipedia.org/wiki/Convolution#Discrete_convolution). That is to say, one way to perform a discrete convolution would be to convert all the terms to coefficients of two polynomials, multiply them together and distribute, and then convert back to their original form. In this sense, everything that can be modelled as a discrete convolution can be modelled using these polynomials and convoluted using the FFT.

In image processing, anything that needs a [kernel](https://en.wikipedia.org/wiki/Kernel_(image_processing)) can be calculated using the FFT.

# The results
So to count the ways to give change, I create my generating functions, then I multiply them all together using the FFT, and then I plug the results into that function I gave above to get the answer.

### Does it work?
Sometimes. I hit a very unexpected roadblock for the FFT implementation. Well before I begin to lose precision from floating point operations, the coefficients begin to lose precision from floating point rounding error, even when using doubles. This will happen once you get a number [too large to store in the 52 bit mantissa of a double without rounding](https://stackoverflow.com/questions/3793838/which-is-the-first-integer-that-an-ieee-754-float-is-incapable-of-representing-e), which happens to be around 10 quadrillion for doubles. The coefficients of these polynomials were getting *that big* for certain currencies.

What's neat is that there's a lot of currencies that do still work with it, like what I call the circulation US dollar, which is all denominations minus the 50 cent and $2 bill.

### Is it faster?
The FFT is generally faster than the dynamic programming approach, but perhaps disappointingly, it's not faster than an optimized "dumb multiply" which can take advantage of the regular structure and predictable intervals of the polys this uses.

Theoretically, as the multiplies get bigger, and the number of coefficients goes up, the FFT should win out in the end. Ultimately, I've found the coefficients get too large to represent too quickly before this can have an impact.

### Is it more space efficient than the dynamic programming?
This is unfortunately a no. Surprisingly, my first unoptimized attempt at the dynamic programming solution (using an STL map that stores keys and values) has a smaller memory footprint than the space it takes to keep one polynomial. The relative sizes can vary (thanks to that GCD optimization I discussed above), but in general it's smaller or very close.

For example, the "circulation US dollar" has 16,272 coefficients that use ~130 KB, while the euro with its mighty 2 cent coin invalidates the GCD optimization and as a result has over 1.4 million coefficients that use over 11 MB of space. This utilization is independent of the actual amount to give change for.

Meanwhile the Dynamic programming utilizations which are dependent on the input change amount is in the ballpark of 250 KB for that USD case, and less than 500 KB for the euro when pushed to the upper end of what unsigned 64 bit integers can reasonably handle for the output count.

The polynomial *could* be optimized to be half as small, (because all input polynomials are symmetric, the output polynomial is symmetric too), but the DP version is just so much smaller than it and can itself be optimized using arrays to make it much smaller in the long run.

# Conclusion
This little experiment was always going to be more about the journey than the destination. The "making change" is a toy problem meant to teach dynamic programming, and this was only a juvenile attempt to solve it the 'wrong' way.

Along the way I followed in the footsteps of mathematicians and computer scientists before me, and rediscovered many of the uses and use cases for techniques and algorithms I once thought monolithic and intimidating. What I learned from this is much greater than ways to give change, and I hope to pass it on to others.

# lda1pfs
An implementation of latent Dirichlet allocation using sequential Bayesian inference

Formerly hosted at code.google.com/p/lda1pfs

This applies the inference style presented in 1PFS (1 Pass Filtering with Shrinkage), Suhrid Balakrishnan and David Madigan, http://paul.rutgers.edu/~suhrid/code/1PFS/ to LDA (latent Dirichlet allocation, http://en.wikipedia.org/wiki/Latent_Dirichlet_allocation ), a graphical model typically used for topic discovery in text.

This implementation uses a Gibbs sampler for the initial phase of 1PFS, where a likely parameter sample set is generated from a fixed-size subset of the entire data set (which may be presented as a stream without a priori bound) and then maintains that sample set in the final phase with a particle filter with shrinkage. For the initial phase, I closely followed existing Gibbs samplers for LDA, which are quite common, and for the final phase, I closely followed the method in 1PFS. To do this, the theta matrix parameter (representing topic distributions within documents, for text) can no longer be stored, so it is re-estimated using phi (topic distributions for words) and the documents as needed using a quick miniature round of Gibbs sampling.

The main advantage of this style of inference is that there is no longer any direct dependence on the length of the data set (number of documents, in the case of text) but there is still indirect dependence through vocabulary size and perhaps if adaptive phi sample sizes were used, also there. However in either of these cases the dependence should be much less than the linear dependence in other inference schemes of which I'm aware.

As far as I can tell, the defining traits of the model survived this intact, but I haven't thought too carefully about this and I could always be wrong. Practically speaking, I would guess that one should expect good results only if the initial Gibbs sampler is able to take into account much of the vocabulary space, since everything depends so heavily on getting a reasonably comprehensive phi estimate before the particle filtering begins, and one will typically be tempted to choose a smallish sample size for the particle filter when not running the samples' computations in parallel across many machines, or maybe even then, in order to be able to stream new data in quickly. In the case of text, where vocabulary size going as the square root of corpus size is the rule of thumb I use, it should usually be the case that much of the vocabulary will come into evidence early on, and it was in my very limited experiments.

It seems to work pretty well at a glance. I wrote it in the course of collaborative filtering work and tested it on some internal data right after an apparently successful run against 20-newsgroups ( http://people.csail.mit.edu/jrennie/20Newsgroups/ ) but for that purpose I quickly replaced it with a more ad-hoc sparse-sampling-on-a-graph algorithm that scaled better yet (that is, well enough, as opposed to not quite well enough, all things considered), and in my own time I've moved on to implementing some older ideas of mine from spring 2008 about bias-optimal Bayesian sequence prediction and planning with respect to the universal distribution, and don't personally plan to invest too much more in this.

Nonetheless, I'm hoping it will be useful or instructive to someone, and if you plan to use it or would like help using it, please let me know. The main to-do item I have for this is to parallelize the particle filter, which I'd be happy to discuss, though I don't see myself getting around to it any time soon. I might write a demonstration script in Scala sometime, for fun.

Oh, and the documentation. There is some in the doc directory. Just a little. To get you started.

This project is provided here under the GNU General Public License v3, except where otherwise noted (the Mersenne Twister bits come to mind, which I might eliminate in favor of something simpler at some point). 

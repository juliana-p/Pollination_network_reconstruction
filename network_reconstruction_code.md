Network reconstruction
================

<b> Ref: Pereira, Battiston, Ribeiro & Jordán (in prep.) <br> Contact:
Juliana Pereira (<julianapereira.mailto@gmail.com>) <br> Date:
06.Jan.2021 </b> <br> <br>

We show how to use network reconstruction based on stochastic block
models on a bipartite network of species interactions (here,
plant-pollinator). <br> The network reconstruction is done in python. I
have used R for the preparation of the data and checking of the results.
If you are not familiar with python tools, the python section can be
done in RStudio, using python chunks, although that is slower than
e.g. jupyter notebook.

## Prepare input graph in R

In this section, we will take a raw, incomplete network of species
interactions, and use network inference based on stochastic block models
to reconstruct it.

First, we will prepare the graph in R. Then the network reconstruction
is done in Python, with the package graph-tool. Afterwards, we come back
to R to have a look at the result.

We start by creating a graph from the interaction data we have, as a
binary network (i.e. not accounting for the repetitions of each pair).

``` r
head(data)
```

    ##   from      to       
    ## 1 "PasEdul" "ApiMell"
    ## 2 "ApiMell" "PasMori"
    ## 3 "BacSagi" "XylThye"
    ## 4 "XylThye" "VerDisc"
    ## 5 "XylThye" "DalFrut"
    ## 6 "XylThye" "DesUnci"

We will also need a data frame of the nodes, informing their mode (who
is plant and who is pollinator):

``` r
head(nodes)
```

    ##            name  mode
    ## AdeBras AdeBras plant
    ## AdeCori AdeCori plant
    ## AdeVisc AdeVisc plant
    ## AecOrna AecOrna plant
    ## AesPauc AesPauc plant
    ## AgeCony AgeCony plant

Let’s build the graph and add the mode as node attribute:

``` r
library(igraph)
g = graph_from_edgelist(data, directed = FALSE)
V(g)$mode = nodes[V(g)$name,2]
g
```

    ## IGRAPH 6043e44 UN-- 1273 2895 -- 
    ## + attr: name (v/c), mode (v/c)
    ## + edges from 6043e44 (vertex names):
    ##  [1] PasEdul--ApiMell ApiMell--PasMori BacSagi--XylThye XylThye--VerDisc
    ##  [5] XylThye--DalFrut XylThye--DesUnci XylThye--IngSess XylThye--PsiBras
    ##  [9] XylThye--CeiSpec XylThye--PseGran XylThye--EucMolu XylThye--GenInfu
    ## [13] XylThye--GueVibu ManFuni--XylChir DesUnci--XylChir IngSess--XylChir
    ## [17] CeiSpec--XylChir EucMolu--XylChir PseGran--XylChir GenInfu--XylChir
    ## [21] XylChir--ChiAlba GueVibu--XylChir XylChir--RanArma CeiSpec--XylIsao
    ## [25] VerDisc--XylPist BacSagi--XylXylo IngSess--XylXylo PsiBras--XylXylo
    ## [29] CeiSpec--XylXylo EucMolu--XylXylo XylXylo--BatAust GenInfu--XylXylo
    ## + ... omitted several edges

Next, let’s format it as a bipartite network. For this, we simply add a
binary node attribute called “type”.

``` r
V(g)$type = 1 * (V(g)$mode == "pollinator")
head(as_data_frame(g, "vertices"))
```

    ##            name       mode type
    ## PasEdul PasEdul      plant    0
    ## ApiMell ApiMell pollinator    1
    ## PasMori PasMori      plant    0
    ## BacSagi BacSagi      plant    0
    ## XylThye XylThye pollinator    1
    ## VerDisc VerDisc      plant    0

We also add an edge attribute `q`, between 0 and 1, which is the
uncertainty of each edge we see in the data. We have tested different
values for our data, and chosen 0.98, because it gave the best result in
the reconstruction (no forbidden edges created, and very few of the
original edges lost).

``` r
E(g)$q = 0.98
head(as_data_frame(g))
```

    ##      from      to    q
    ## 1 PasEdul ApiMell 0.98
    ## 2 ApiMell PasMori 0.98
    ## 3 BacSagi XylThye 0.98
    ## 4 XylThye VerDisc 0.98
    ## 5 XylThye DalFrut 0.98
    ## 6 XylThye DesUnci 0.98

Since our network is a bipartite one, there are forbidden edges,
i.e. edges that we do not want to show up in the reconstruction. To
help prevent them, we want to set the uncertainty about such edges to
zero. For that, we need to include in our graph g all of the possible
forbidden edges, and set their attribute `q` to zero.

For that, we create two complete graphs (i.e. graph where all nodes are
linked directly to all others), one for the plants and one for the
pollinators. We have 584 plants and 689 pollinators in our dataset. So:

``` r
plants_g = make_full_graph(584)
V(plants_g)$name = subset(as_data_frame(g, "vertices"), mode == "plant")$name

pols_g = make_full_graph(689)
V(pols_g)$name = subset(as_data_frame(g, "vertices"), mode == "pollinator")$name
```

Then we create a large graph uniting all three, and set the `q`
attribute to 0 for the forbidden edges:

``` r
G = union(g, plants_g, pols_g)
G
```

    ## IGRAPH fd6e822 UN-B 1273 410147 -- 
    ## + attr: name_2 (g/c), name_3 (g/c), loops_2 (g/l), loops_3 (g/l), mode
    ## | (v/c), type (v/n), name (v/c), q (e/n)
    ## + edges from fd6e822 (vertex names):
    ##  [1] OurHexa--TitDive SolChil--TitDive SolChil--OurHexa SenCras--TitDive
    ##  [5] SenCras--OurHexa SenCras--SolChil GalParv--TitDive GalParv--OurHexa
    ##  [9] GalParv--SolChil GalParv--SenCras EreHier--TitDive EreHier--OurHexa
    ## [13] EreHier--SolChil EreHier--SenCras EreHier--GalParv EmiFosb--TitDive
    ## [17] EmiFosb--OurHexa EmiFosb--SolChil EmiFosb--SenCras EmiFosb--GalParv
    ## [21] EmiFosb--EreHier CoeAcul--CoeSimi AcaPrin--CoeSimi AcaPrin--CoeAcul
    ## [25] CerDarw--CoeSimi CerDarw--CoeAcul CerDarw--AcaPrin ColPetr--CoeSimi
    ## + ... omitted several edges

The forbidden edges will have `q` attribute = NA. Change it to zero:

``` r
E(G)$q[which(is.na(E(G)$q))] = 0
table(as_data_frame(G)$q)
```

    ## 
    ##      0   0.98 
    ## 407252   2895

Now we save the graph in graphml format:

``` r
write_graph(G, "G.graphml", format="graphml")
```

And move to python.

## Network reconstruction in python

We will use the graph-tool package. For instruction on installation and
package documentation, visit: <https://graph-tool.skewed.de/>

Load libraries:

``` python
from graph_tool.all import *
import matplotlib.pyplot as plt
import numpy as np
```

Import the network data:

``` python
G = load_graph("G.graphml")
```

Collect the number of vertices, the number of edges (without taking into
account the forbidden ones\!), the uncertainties for the edges `q`
(saved as the `q` attribute in the graph), and the uncertainties for the
non-edges `q_default` (i.e. all plant-pollinator pairs for which no
interaction was observed).

``` python
N = G.num_vertices()
E = 2895 

q = G.ep.q  # as defined above, 0.98 for observed and 0 for forbidden edges

# The q_default value is chosen to preserve the expected density of the original network
# (see https://graph-tool.skewed.de/static/doc/demos/inference/inference.html#network-reconstruction -> 
# Extraneous error estimates)
# For a bipartite network, use the formula below, where 584 is our number of plants and 689 of pollinator spp
q_default = (E - q.a.sum()) / ((584 * 689) - E)
```

Now the network reconstruction. For details about the functions, see the
graph-tool package documentation. The chunk below may take a few minutes
or longer, depending on the size of the network.

``` python
state = UncertainBlockState(G, q=q, q_default=q_default)
mcmc_equilibrate(state, wait=2000, mcmc_args=dict(niter=10))       
```

The chunk below may take a few hours.

``` python
u = None              # marginal posterior edge probabilities
bs = []               # partitions

def collect_marginals(s):
    global u, bs, cs
    u = s.collect_marginal(u)
    bstate = s.get_block_state()
    bs.append(bstate.levels[0].b.a.copy())

mcmc_equilibrate(state, force_niter=10000, mcmc_args=dict(niter=10),
                    callback=collect_marginals)

eprob = u.ep.eprob
```

`u` is the reconstructed network.

Let’s export the reconstructed network in graphml format:

``` python
u.save("reconstructed.graphml")
```

And go back to R.

## Checking reconstructed network in R

``` r
u = read.graph("reconstructed.graphml",format="graphml")
```

We can copy the attributes of the raw network over to the reconstructed
one (the node id’s should be in the same order, but it is good to check
first):

``` r
# checking the the order of nodes in the 2 graphs is identical
all(as_data_frame(g, what = "vertices")$id == as_data_frame(u, what="vertices")$id) # should be TRUE
```

    ## [1] TRUE

``` r
V(u)$name = V(g)$name # copy species names
V(u)$mode = V(g)$mode # copy mode (plant or pollinator)

head(as_data_frame(u))
```

    ##      from      to count      eprob id
    ## 1 PasEdul ApiMell  9999 1.00000000 e0
    ## 2 PasEdul XylFron  9999 1.00000000 e1
    ## 3 PasEdul TetAngu  9977 0.99779978 e2
    ## 4 PasEdul EugIgni     4 0.00040004 e3
    ## 5 PasEdul VanBraz     9 0.00090009 e4
    ## 6 PasEdul TriSpin   104 0.01040104 e5

The ‘eprob’ attribute is the posterior probability of the edge according
to the reconstruction. That is, the probability that the link exists. We
have interpreted this value as the probability or potential of
interaction between species.

We can check how many links of the raw network were kept/lost, how many
new ones were created, and how many forbidden links (plant-plant or
pollinator-pollinator) were created:

``` r
count_impossible_links = function(g){
  link_ends = ends(g, E(g))
  mode_from = vertex_attr(g, name = "mode", index = as.character(link_ends[,1]))
  mode_to = vertex_attr(g, name = "mode", index = as.character(link_ends[,2]))
  impossible = length(which(mode_from == mode_to))
  possible = length(which(mode_from != mode_to))
  kept = length(E(intersection(g, gpol)))
  lost = length(E(difference(gpol,g)))
  created = length(E(difference(g, gpol)))
  return(c("impossible"=impossible,"possible"=possible,"kept"=kept,"lost"=lost,"created"=created))
}

count_impossible_links(u)
```

    ## impossible   possible       kept       lost    created 
    ##          0      49406       2893          2      46513

The ideal is to have zero forbidden links and zero or very few lost
links. If that doesn’t happen, try different `q` values for the observed
edges.

This then is our network of potential of species interactions.
Multiplying it for each local network of co-occurrence probabilities, we
obtained our local pollination networks.

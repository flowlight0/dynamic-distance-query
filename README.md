Fully Dynamic Distance Query
==============================
This is an implementation of a distance query algorithm on massive dynamic networks (CIKM '16). This implementation enables you to efficiently compute exact point-to-point distance on dynamic large complex networks (social networks, web graphs, co-author networks). If you have interest in our method and experimental results, please see [our paper](http://dl.acm.org/citation.cfm?doid=2983323.2983731).

## Usage
### CUI
We can use our CUI as the following: 

    $ ./waf configure
    $ ./waf
    $ ./bin/dynamic_distance_query --graph_file sample/graph.txt --query_file sample/query.txt --undirected
    5
    2
    7
    1
    2
    6
    Invalid op DUMMY, skipping it
You can see the format of input graph file and query file from https://github.com/flowlight0/dynamic-distance-query/blob/master/sample/graph.txt and https://github.com/flowlight0/dynamic-distance-query/blob/master/sample/query.txt, respectively. 

### Benchmarker
Theare are two benchmarkers for examining the performance of our algorithm. https://github.com/flowlight0/dynamic-distance-query/blob/master/src/cui/benchmark_various_search_strategy_main.cpp is a benchmarker for comparing the querying efficiency of several search strategies in static settings. https://github.com/flowlight0/dynamic-distance-query/blob/master/src/cui/benchmark_dynamic_update_main.cpp is a benchmarker for the performance of our method with various index size reduction parameter in dynamic settings. 
For example, we can run `./bin/benchmark_various_search_strategy` and `./bin/benchmark_dynamic_update` with the following commands. 

    $ ./bin/benchmark_various_search_strategy --graph_file=data/samplegraph.txt --num_bfs_queries=100 --num_bibfs_queries=1000 --undirected
    $ ./bin/benchmark_dynamic_update --graph_file=data/samplegraph.txt --undirected --num_updates=5


These benchmark scripts generate json files containing results of benchmarks. Since two benchmarkers are only assuming undirected graphs. If you forget to add `--undirected` option, they might give unexpected results. 

## Reference 
Takanori Hayashi, Takuya Akiba, and Ken-ichi Kawarabayashi. [**Fully Dynamic Shortest-Path Distance Query Acceleration on Massive Networks**](http://dl.acm.org/citation.cfm?doid=2983323.2983731).  [*CIKM'16*](http://www.cikm2016.org/)

## Contact 
If you have questions or find errors, please contact to me (flowlight0 at gmail.com)

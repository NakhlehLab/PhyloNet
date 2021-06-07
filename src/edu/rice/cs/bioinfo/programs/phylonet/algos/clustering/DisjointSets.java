package edu.rice.cs.bioinfo.programs.phylonet.algos.clustering;

import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: zhujiafan
 * Date: 10/8/16
 * Time: 11:37 AM
 * To change this template use File | Settings | File Templates.
 */
public class DisjointSets<T> {
    int _n;
    List<T> _elements;
    Map<T, Integer> _indices;
    int[] _parents;
    int[] _sizes;

    public DisjointSets(Set<T> elements) {
        _n = elements.size();
        _parents = new int[_n];
        _sizes = new int[_n];
        _elements = new ArrayList<>();
        _indices = new HashMap<>();
        for(T element : elements) {
            _indices.put(element, _elements.size());
            _elements.add(element);
        }
        for(int i = 0 ; i < _n ; i++) {
            _parents[i] = i;
            _sizes[i] = 1;
        }
    }

    private int find(int i) {
        while(i != _parents[i]) {
            _parents[i] = _parents[_parents[i]];
            i = _parents[i];
        }
        return i;
    }

    private void union(int p, int q) {
        int i = find(p);
        int j = find(q);
        if(i == j) return;
        if(_sizes[i] < _sizes[j]) {
            _parents[i] = j;
            _sizes[j] += _sizes[i];
        } else {
            _parents[j] = i;
            _sizes[i] += _sizes[j];
        }
    }

    public void union(T a, T b) {
        int p = _indices.get(a);
        int q = _indices.get(b);
        union(p, q);
    }

    public void union(Set<T> elements) {
        if(elements.size() <= 1) return;
        T prev = null;
        for(T e : elements) {
            if(prev != null) {
                union(prev, e);
            }
            prev = e;
        }
    }

    public Set<Set<T>> getDisjointSets() {
        Map<Integer, Set<T>> indexToSets = new HashMap<>();
        for(T e : _indices.keySet()) {
            int i = _indices.get(e);
            int p = find(i);
            if(!indexToSets.containsKey(p)) {
                indexToSets.put(p, new HashSet<T>());
            }
            indexToSets.get(p).add(e);
        }
        return new HashSet<>(indexToSets.values());
    }

}

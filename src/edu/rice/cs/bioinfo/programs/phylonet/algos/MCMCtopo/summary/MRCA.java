package edu.rice.cs.bioinfo.programs.phylonet.algos.MCMCtopo.summary;

import edu.rice.cs.bioinfo.library.programming.extensions.java.lang.iterable.IterableHelp;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.Networks;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;

/**
 * To summarize MRCA (most recent common ancestor) clusters given samples
 *
 * Created by wendingqiao on 3/6/15.
 */
public class MRCA {

    private Map<Network, Double> clusters;
    private double total = 0.0;

    public MRCA(String file) {
        this.clusters = new HashMap<>();
        mrcaSumFile(file);
    }

    public void report() {
        for(Network key : clusters.keySet()) {
            double value = clusters.get(key) / total;
            System.out.println(value + "  " + key.toString());
        }
    }

    private void mrcaSumFile(String file) {
        try{
            BufferedReader in = new BufferedReader(new FileReader(file));
            String s;
            while((s = in.readLine()) != null) {
                String[] ss = s.split("\\s+");
                if(ss.length == 5) {
                    int size = Integer.parseInt(ss[0].substring(0, ss[0].length() - 1));
                    double averagePost = Double.parseDouble(ss[2].substring(0, ss[2].length() - 1));
                    double tmp = averagePost * size;
                    Network net = Networks.readNetwork(ss[3]);

                    List<Network> previousClusters = new ArrayList<>();
                    for(Object n : net.getNetworkNodes()) {
                        NetNode node = (NetNode) n;
                        Network cluster = getCluster(node);
                        boolean found = false;
                        for(Network pc : previousClusters) {
                            if(Networks.hasTheSameTopology(pc, cluster)) {
                                found = true;
                                break;
                            }
                        }
                        if(found) continue;
                        previousClusters.add(cluster);
                        for(Network key : clusters.keySet()) {
                            if(Networks.hasTheSameTopology(cluster, key)) {
                                clusters.put(key, clusters.get(key) + tmp);
                                found = true;
                                break;
                            }
                        }
                        if(!found) clusters.put(cluster, tmp);
                    }
                    total += tmp;
                }
            }
            in.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private Network getCluster(NetNode node) {
        Map<NetNode, BniNetNode> map = new HashMap<>();
        while(true) {
            List<NetNode> children = IterableHelp.toList(node.getChildren());
            if(children != null && children.size() == 1) {
                node = children.get(0);
            } else break;
        }
        Queue<NetNode> queue = new LinkedList<>();
        queue.add(node);
        while(!queue.isEmpty()) {
            NetNode cur = queue.poll();
            BniNetNode curNew = new BniNetNode();
            curNew.setName(cur.getName());
            map.put(cur, curNew);
            for(Object child : cur.getChildren())  {
                queue.add((NetNode) child);
            }
        }
        Network cluster = new BniNetwork(map.get(node));
        queue.add(node);
        while(!queue.isEmpty()) {
            NetNode cur = queue.poll();
            for(Object p : cur.getParents())  {
                NetNode par = (NetNode) p;
                if(map.containsKey(par)) {
                    map.get(par).adoptChild(map.get(cur), 1);
                }
            }
            for(Object child : cur.getChildren())  {
                queue.add((NetNode) child);
            }
        }
        for(Object n : cluster.dfs()) {
            NetNode cur = (NetNode) n;
            if(cur.isRoot() || cur.isNetworkNode()) continue;
            List<NetNode> parents = IterableHelp.toList(cur.getParents());
            List<NetNode> children = IterableHelp.toList(cur.getChildren());
            if(parents.size() == 1 && children.size() == 1) {
                NetNode parent = parents.get(0);
                for(NetNode c : children) {
                    cur.removeChild(c);
                    parent.adoptChild(c, 1);
                }
                parent.removeChild(cur);
            }
        }
        return cluster;
    }
}

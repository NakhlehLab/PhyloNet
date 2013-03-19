/*
 * Copyright (c) 2012 Rice University.
 *
 * This file is part of PhyloNet.
 *
 * PhyloNet is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PhyloNet is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PhyloNet.  If not, see <http://www.gnu.org/licenses/>.
 */

package edu.rice.cs.bioinfo.programs.phylonet.algos.simulator;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: Matt
 * Date: 11/2/11
 * Time: 10:47 AM
 * To change this template use File | Settings | File Templates.
 */
public class SimGTInNetwork2010
{
    private double t1, t2, gamma;
	private int n;
	private String output;

	public double g22(double t){
		return Math.exp((-1)*t);
	}

	public double g21(double t){
		return 1-Math.exp((-1)*t);
	}

	public double g32(double t){
		return 3.0/2*Math.exp((-1)*t) - 3.0/2*Math.exp((-3)*t);
	}

	public double g33(double t){
		return Math.exp((-3)*t);
	}

	public double g31(double t){
		return 1-g32(t)-g33(t);
	}




	public double exprnd(double mean){
		 double u = Math.random();
	     return -(mean * Math.log(u));
	}


	public String[] generateGTs(double t1, double t2, double gamma, int n){
		double x1,x2,x3;
		String gtexp = "";
		String[] gtlist = new String[n];

		for(int i=0; i<n; i++){
		    x1 = exprnd(1);
		    if(x1 <= t2){  //B and C coalesce more recently than hybridization
		        double u = Math.random();  //u decides whether (B,C) goes left or right
		        x2 = exprnd(1);
		        if(x2 <= t1){  //(B,C) and A(or D) coalesce more recently than root
		            x3 = exprnd(1);
		            if(u <= gamma){   //(B,C) goes left
		                gtexp = "(((B:"+ x1 + ",C:" + x1 + "):" + (t2-x1+x2) + ",A:" +(t2+x2) + "):" + (t1-x2+x3) + ",D:" + (t1+t2+x3) + ");";
		                //senario = 'B and C coalesce within t2, then go left and coalesce with A within t1';
		            }
		            else{
		            	gtexp = "(((B:"+ x1 + ",C:" + x1 + "):" + (t2-x1+x2) + ",D:" +(t2+x2) + "):" + (t1-x2+x3) + ",A:" + (t1+t2+x3) + ");";
		                //senario = 'B and C coalesce within t2, then go right and coalesce with D within t1';
		            }
		        }else{
		            double p32 = Math.random();
		            x2 = exprnd(1.0/3);
		            x3 = exprnd(1);
		            if(p32 <= 1.0/3){  //(B,C) and A coalesce first above root
		                gtexp = "(((B:" + x1 + ",C:" + x1 + "):" + (t2-x1+t1+x2) + ",A:" + (t2+t1+x2) + "):" + x3 + ",D:" + (t2+t1+x2+x3) + ");";
		                //senario = 'B and C coalesce within t2, then (B,C) coalesce with A above the root, then with D';
		            }
		            else if(p32 <= 2.0/3){   //(B,C) and D coalesce first above root
		            	gtexp = "(((B:" + x1 + ",C:" + x1 + "):" + (t2-x1+t1+x2) + ",D:" + (t2+t1+x2) + "):" + x3 + ",A:" + (t2+t1+x2+x3) + ");";
		                //senario = 'B and C coalesce within t2, then (B,C) coalesce with D above the root, then with A';
		            }
		            else{  //A and D coalesce first above root
		                gtexp = "((B:" + x1 + ",C:" + x1 + "):" + (t2-x1+t1+x2+x3) + ",(A:" + (t2+t1+x2) + ",D:" + (t2+t1+x2) + "):" + x3 + ");";
		                //senario = "B and C coalesce within t2, then coalesce with (A,D) above the root";
		            }
		        }
		    }
		    else{
		        double u1 = Math.random();  //u1 decides whether B goes left or right
		        double u2 = Math.random();  //u2 decides whether C goes left or right
		        boolean allaboveroot = false;
		        if((u1 <= gamma) && (u2 <= gamma)){  //both B and C goes left
		            x1 = exprnd(1.0/3);
		            if(x1 <= t1){
		                double p32 = Math.random();
		                String rest;
		                if(p32 <= 1.0/3){  //B and C coalesce first
		                    gtexp = "(B:" + (t2+x1) + ",C:" + (t2+x1) + ")";
		                    rest = "A:";
		                }
		                else if(p32 <= 2.0/3){  //A and C coalesce first
		                    gtexp = "(A:" + (t2+x1) + ",C:" + (t2+x1) + ")";
		                    rest = "B:";
		                }
		                else{  //A and B coalesce first
		                    gtexp = "(A:" + (t2+x1) + ",B:" + (t2+x1) + ")";
		                    rest = "C:";
		                }
		                x2 = exprnd(1);
		                if(x2 <= t1-x1){  //A, B and C all coalesce below the root
		                    x3 = exprnd(1);
		                    gtexp = "((" + gtexp + ":" + x2 + "," + rest + (t2+x1+x2) + "):" + (t1-x1-x2+x3) + ",D:"+ (t1+t2+x3) + ");";
		                    //senario = "both B and C goes left, A, B and C all coalesce below the root";
		                }
		                else{
		                    p32 = Math.random();
		                    x2 = exprnd(1.0/3);
		                    x3 = exprnd(1);
		                    if(p32 <= 1.0/3)  //A, B and C coalesce first above the root
		                        gtexp = "((" + gtexp + ":" + (t1-x1+x2) + "," + rest + (t2+t1+x2) + "):" + x3 + ",D:" + (t1+t2+x2+x3) + ");";
		                    else if(p32 <= 2.0/3)  //the coalesced node and D coalesce first above the root
		                        gtexp = "((" + gtexp + ":" + (t1-x1+x2) + ",D:" + (t2+t1+x2) + "):" + x3 + "," + rest + (t1+t2+x2+x3) + ");";
		                    else
		                        gtexp = "(" + gtexp + ":" + (t1-x1+x2+x3) + ",(D:" + (t2+t1+x2) + "," + rest + (t2+t1+x2) + "):" + x3 + ");";
		                }
		                //senario = "both B and C goes left, two of A, B and C coalesce below the root";
		            }
		            else{
		                allaboveroot = true;
		            }
		        }
		        else if((u1 > gamma) && (u2 > gamma)){  //both B and C goes right
		            x1 = exprnd(1.0/3);
		            String rest;
		            if(x1 <= t1){
		                double p32 = Math.random();
		                if(p32 <= 1.0/3){  //B and C coalesce first
		                    gtexp = "(B:" + (t2+x1) + ",C:" + (t2+x1) + ")";
		                    rest = "D:";
		                }
		                else if(p32 <= 2.0/3){  //D and C coalesce first
		                    gtexp = "(D:" + (t2+x1) + ",C:" + (t2+x1) + ")";
		                    rest = "B:";
		                }
		                else{  //%D and B coalesce first
		                    gtexp = "(D:" + (t2+x1) + ",B:" + (t2+x1) + ")";
		                    rest = "C:";
		                }
		                x2 = exprnd(1);
		                if(x2 <= t1-x1){  //D, B and C all coalesce below the root
		                    x3 = exprnd(1);
		                    gtexp = "((" + gtexp + ":" + x2 + "," + rest + (t2+x1+x2) + "):" + (t1-x1-x2+x3) + ",A:" + (t1+t2+x3) + ");";
		                    //senario = "both B and C goes right, D, B and C all coalesce below the root";
		                }
		                else{
		                    p32 = Math.random();
		                    x2 = exprnd(1.0/3);
		                    x3 = exprnd(1);
		                    if(p32 <= 1.0/3)  //A, B and C coalesce first above the root
		                        gtexp = "((" + gtexp + ":" + (t1-x1+x2) + "," + rest + (t2+t1+x2) + "):" + x3 + ",A:" + (t1+t2+x2+x3) + ");";
		                    else if(p32 <= 2.0/3)   // the coalesced node and A coalesce first above the root
		                        gtexp = "((" + gtexp + ":" + (t1-x1+x2) + ",A:" + (t2+t1+x2) + "):" + x3 + "," + rest + (t1+t2+x2+x3) + ");";
		                    else
		                        gtexp = "(" + gtexp + ":" + (t1-x1+x2+x3) + ",(A:" + (t2+t1+x2) + "," + rest + (t2+t1+x2) + "):" + x3 + ");";
		                    //senario = "both B and C goes right, two of D, B and C coalesce below the root";
		                }
		            }
		            else
		                allaboveroot = true;
		        }
		        else if((u1 <= gamma) && (u2 > gamma)){
		            int coal = 0;
		            x1 = exprnd(1);
		            if(x1<=t1){  //B and A coalesce below the root
		                coal = 1;
		            }
		            x2 = exprnd(1);
		            if(x2<=t1){  //C and D coalesce below the root
		                if(coal ==1)
		                    coal = 3;
		                else
		                    coal = 2;
		            }
		            if(coal==1){
		                x2 = exprnd(1.0/3);
		                x3 = exprnd(1);
		                double p32 = Math.random();
		                if(p32<=1.0/3)
		                    gtexp = "(((A:" + (t2+x1) + ",B:" + (t2+x1) + "):" + (t1-x1+x2) + ",C:" + (t2+t1+x2) + "):" + x3 + ",D:" + (t2+t1+x2+x3) + ");";
		                else if(p32<=2.0/3)
		                    gtexp = "(((A:" + (t2+x1) + ",B:" + (t2+x1) + "):" + (t1-x1+x2) + ",D:" + (t2+t1+x2) + "):" + x3 + ",C:" + (t2+t1+x2+x3) + ");";
		                else
		                    gtexp = "((A:" + (t2+x1) + ",B:" + (t2+x1) + "):" + (t1-x1+x2+x3) + ",(C:" + (t2+t1+x2) + ",D:" + (t2+t1+x2) + "):" + x3 + ");";
		            }
		                //senario = "B goes left and C goes right; B and A coalesce below the root";
		            else if(coal==2){
		                x1 = x2;
		                x2 = exprnd(1.0/3);
		                x3 = exprnd(1);
		                double p32 = Math.random();
		                if(p32<=1.0/3)
		                    gtexp = "(((C:" + (t2+x1) + ",D:" + (t2+x1) + "):" + (t1-x1+x2) + ",A:" + (t2+t1+x2) + "):" + x3 + ",B:" + (t2+t1+x2+x3) + ");";
		                else if(p32<=2.0/3)
		                    gtexp = "(((C:" + (t2+x1) + ",D:" + (t2+x1) + "):" + (t1-x1+x2) + ",B:" + (t2+t1+x2) + "):" + x3 + ",A:" + (t2+t1+x2+x3) + ");";
		                else
		                    gtexp = "((C:" + (t2+x1) + ",D:" + (t2+x1) + "):" + (t1-x1+x2+x3) + ",(A:" + (t2+t1+x2) + ",B:" + (t2+t1+x2) + "):" + x3 + ");";
		            }
		                //senario = "B goes left and C goes right; C and D coalesce below the root";
		            else if(coal==3){
		                x3 = exprnd(1);
		                gtexp = "((A:" + (t2+x1) + ",B:" + (t2+x1) + "):" + (t1-x1+x3) + ",(C:" + (t2+x2) + ",D:" + (t2+x2) + "):" + (t1-x2+x3) + ");";
		                //senario = "B goes left and C goes right; B and A coalesce below the root; C and D coalesce below the root";
		            }
		            else
		                allaboveroot = true;
		        }
		        else{
		            int coal = 0;
		            x1 = exprnd(1);
		            if(x1<=t1){  //C and A coalesce below the root
		                coal = 1;
		            }
		            x2 = exprnd(1);
		            if(x2<=t1){ //B and D coalesce below the root
		                if(coal ==1)
		                    coal = 3;
		                else
		                    coal = 2;
		            }
		            if(coal==1){
		                x2 = exprnd(1.0/3);
		                x3 = exprnd(1);
		                double p32 = Math.random();
		                if(p32<=1.0/3)
		                    gtexp = "(((A:" + (t2+x1) + ",C:" + (t2+x1) + "):" + (t1-x1+x2) + ",B:" + (t2+t1+x2) + "):" + x3 + ",D:" + (t2+t1+x2+x3) + ");";
		                else if(p32<=2.0/3)
		                    gtexp = "(((A:" + (t2+x1) + ",C:" + (t2+x1) + "):" + (t1-x1+x2) + ",D:" + (t2+t1+x2) + "):" + x3 + ",B:" + (t2+t1+x2+x3) + ");";
		                else
		                    gtexp = "((A:" + (t2+x1) + ",C:" + (t2+x1) + "):" + (t1-x1+x2+x3) + ",(B:" + (t2+t1+x2) + ",D:" + (t2+t1+x2) + "):" + x3 + ");";
		                //senario = "B goes right and C goes left; B and D coalesce below the root";
		            }
		            else if(coal==2){
		                x1 = x2;
		                x2 = exprnd(1.0/3);
		                x3 = exprnd(1);
		                double p32 = Math.random();
		                if(p32<=1.0/3)
		                    gtexp = "(((B:" + (t2+x1) + ",D:" + (t2+x1) + "):" + (t1-x1+x2) + ",A:" + (t2+t1+x2) + "):" + (x3) + ",C:" + (t2+t1+x2+x3) + ");";
		                else if(p32<=2.0/3)
		                    gtexp = "(((B:" + (t2+x1) + ",D:" + (t2+x1) + "):" + (t1-x1+x2) + ",C:" + (t2+t1+x2) + "):" + (x3) + ",A:" + (t2+t1+x2+x3) + ");";
		                else
		                    gtexp = "((B:" + (t2+x1) + ",D:" + (t2+x1) + "):" + (t1-x1+x2+x3) + ",(A:" + (t2+t1+x2) + ",C:" + (t2+t1+x2) + "):" + (x3) + ");";
		                //senario = "B goes right and C goes left; B and D coalesce below the root";
		            }
		            else if(coal==3){
		                x3 =exprnd(1);
		                gtexp = "((A:" + (t2+x1) + ",C:" + (t2+x1) + "):" + (t1-x1+x3) + ",(B:" + (t2+x2) + ",D:" + (t2+x2) + "):" + (t1-x2+x3) + ");";
		                //senario = "B goes right and C goes left; C and A coalesce below the root; B and D coalesce below the root";
		            }
		            else
		                allaboveroot = true;
		        }

		        if(allaboveroot){
		            //senario = "all lineages coalesce above the root";
		            x1 = exprnd(1.0/6);
		            x2 = exprnd(1.0/3);
		            x3 = exprnd(1);
		            List<String> nodes = new ArrayList<String>();
		            nodes.add("A:");
		            nodes.add("B:");
		            nodes.add("C:");
		            nodes.add("D:");
		            double p41 = Math.random();
		            String node1;
		            if(p41<=1.0/4){
		            	node1 = nodes.get(0);
		            	nodes.remove(0);
		            }
		            else if(p41<=2.0/4){
		            	node1 = nodes.get(1);
		            	nodes.remove(1);
		            }
		            else if(p41<=3.0/4){
		            	node1 = nodes.get(2);
		            	nodes.remove(2);
		            }
		            else{
		            	node1 = nodes.get(3);
		            	nodes.remove(3);
		            }
		            double p31 = Math.random();
		            String node2;
		            if(p31<=1.0/3){
		            	node2 = nodes.get(0);
		            	nodes.remove(0);
		            }
		            else if(p31<=2.0/3){
		            	node2 = nodes.get(1);
		            	nodes.remove(1);
		            }
		            else{
		            	node2 = nodes.get(2);
		            	nodes.remove(2);
		            }
		            String newnode = "(" + node1 + (t1+t2+x1) + "," + node2 + (t1+t2+x1) + "):";
		            nodes.add(newnode);
		            double p32 = Math.random();
		            if(p32<=1.0/3)
		                gtexp = "((" + nodes.get(0) + x2 + "," + nodes.get(1) + (t1+t2+x1+x2) + "):" + x3 + "," + nodes.get(2) + (t1+t2+x1+x2+x3) + ")";
		            else if(p32<=2.0/3)
		                gtexp = "((" + nodes.get(0) + x2 + "," + nodes.get(2) + (t1+t2+x1+x2) + "):" + x3 + "," + nodes.get(1) + (t1+t2+x1+x2+x3) + ")";
		            else
		            	gtexp = "((" + nodes.get(1) + x2 + "," + nodes.get(2) + (t1+t2+x1+x2) + "):" + x3 + "," + nodes.get(0) + (t1+t2+x1+x2+x3) + ")";
		        }
		    }

		    gtlist[i] = gtexp;
		}
		return gtlist;
	}
}

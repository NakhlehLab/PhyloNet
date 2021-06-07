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

package edu.rice.cs.bioinfo.programs.phylonet.structs.network.io;

import edu.rice.cs.bioinfo.programs.phylonet.structs.network.NetNode;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.Network;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.model.bni.BniNetwork;
import edu.rice.cs.bioinfo.programs.phylonet.structs.network.util.NetNodes;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;
import java.util.LinkedList;
import java.util.List;

/**
 * This class is for reading and constructing a network from a tree and a set of
 * HGT events. A tree is in the Newick format, and each event is of the form
 * source -> destination, where source and destionation are names of the heads of
 * the donor and receipent edges.
 */
public class HgtReader<T> {
	public HgtReader(Reader r)
	{
		_reader = r;
	}

	/**
	 * This function builds a network from a tree and a set of HGT events. The network
	 * read is stored in the data member _net.
	 */
	public Network<T> readHgt()
	{
		// Read input and add a name to the tree to comply with eNewick format.
		getInput();

		// Parse the main tree.
		ExNewickReader<T> enr = new ExNewickReader<T>(new StringReader(_first_line));

        try
        {
            _net = enr.readNetwork();
        }
        catch(Exception e)
        {
            throw new RuntimeException(e);
        }

		// Add events to the network.
		for (String str : _next_lines) {
			String names[] = parseEvent(str);
			NetNode<T> source = _net.findNode(names[0]);
			NetNode<T> dest = _net.findNode(names[1]);

			if (source == null || dest == null) {
				throw new RuntimeException("Event " + names[0] + " -> " + names[1] + " cannot be added to the network");
			}
			else {
				if (!source.isRoot()) {
					source = NetNodes.breakEdge(source.getParents().iterator().next(), source, NetNode.NO_DISTANCE);
					dest = NetNodes.breakEdge(dest.getParents().iterator().next(), dest, NetNode.NO_DISTANCE);
					source.adoptChild(dest, NetNode.NO_DISTANCE);
				}
				else {
					Network<T> newNet = new BniNetwork<T>();
					newNet.createRoot("ROOT");
					newNet.getRoot().adoptChild(source, NetNode.NO_DISTANCE);
					dest = NetNodes.breakEdge(dest.getParents().iterator().next(), dest, NetNode.NO_DISTANCE);
					source.adoptChild(dest, NetNode.NO_DISTANCE);

					_net = newNet;
				}
			}
		}

		return _net;
	}

	/**
	 * This function returns the name of the source and destination node of an event.
	 *
	 * @param line: An event in the form "source -> destionation"
	 * @return: An array of names, with the names of the source and destination in
	 * indices 0 and 1.
	 */
	private String[] parseEvent(String line)
	{
		String DELIM = "->";
		int pos = line.indexOf(DELIM);

		if (pos < 0) {
			throw new RuntimeException("Cannot parse line " + line + " into an event");
		}
		else {
			String names[] = new String[2];
			names[0] = line.substring(0, pos).trim();
			names[1] = line.substring(pos + DELIM.length(), line.length()).trim();

			return names;
		}
	}

	/**
	 * This function reads the input from the reader. Input strings are stored in
	 * _first_line and _next_lines.
	 */
	private void getInput()
	{
		BufferedReader br = new BufferedReader(_reader);

		try {
			// Read the first non-blank line as the main tree.
			while ((_first_line = br.readLine()) != null) {
				_first_line = _first_line.trim();
				if (_first_line.length() != 0) {
					break;
				}
			}

			// Read subsequent non-blank lines as subtrees.
			String str;
			_next_lines = new LinkedList<String>();

			while ((str = br.readLine()) != null) {
				str = str.trim();
				if (str.length() != 0) {
					_next_lines.add(str);
				}
			}
		}
		catch (IOException e) {
			System.err.print(e.getMessage());
			throw new RuntimeException("Cannot read the network");
		}
	}

	// Data members.
	private String _first_line;			// Stores the main tree.
	private List<String> _next_lines;	// Stores events.
	private Network<T> _net;	// The network to be constructed.
	private Reader _reader;		// Input stream.
}

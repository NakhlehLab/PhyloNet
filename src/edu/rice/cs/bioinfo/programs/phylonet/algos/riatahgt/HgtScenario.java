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

package edu.rice.cs.bioinfo.programs.phylonet.algos.riatahgt;

import java.io.Serializable;
import java.util.LinkedList;
import java.util.List;

public class HgtScenario implements Serializable {
	/* Methods of ExHgtScenario */
	public HgtScenario()
	{
		_events = new LinkedList<HgtEvent>();
	}

	/**
	 * The method prints events in the list of events. Bad events are printed in
	 * the form of refinment on nodes. Violated events have an additional explanation
	 * [time violation?].
	 */
	public String toString()
	{
		String str = new String();
		boolean hasBad = false;

		for (HgtEvent event : _events) {
			if (!event.isBad()) {
				str += event.toString();
				if (event.isViolated()) {
					str += " [time violation?]\n";
				}
				else {
					str += "\n";
				}
			}
			else {
				hasBad = true;
			}
		}

		if (hasBad) {
			str += "[refine nodes: ";
			for (HgtEvent event : _events) {
				if (event.isBad()) {
					str += event.getSourceEdge().getName() + ", ";
				}
			}

			str = str.substring(0, str.length() - 2) + "]\n";
		}

		return str;
	}

	public String toString(int tab_level)
	{
		String str = new String();
		boolean hasBad = false;

		for (HgtEvent event : _events) {
			if (!event.isBad()) {
				str += printTabs(tab_level) + event.toString();
				if (event.isViolated()) {
					str += " [time violation?]\n";
				}
				else {
					str += "\n";
				}
			}
			else {
				hasBad = true;
			}
		}

		if (hasBad) {
			str += printTabs(tab_level) + "[refine nodes: ";
			for (HgtEvent event : _events) {
				if (event.isBad()) {
					str += event.getSourceEdge().getName() + ", ";
				}
			}

			str = str.substring(0, str.length() - 2) + "]\n";
		}

		return str;
	}

	private String printTabs(int tab_level)
	{
		String tabs = "";

		for (int i = 0; i < tab_level; i++) {
			tabs += "\t";
		}

		return tabs;
	}

	/**
	 * This function counts the number of bad events. An event x -> y is bad
	 * if x is the parent of y.
	 *
	 * @return: The number of bad events.
	 */
	public int countBadEvents()
	{
		if (_events == null) {
			return 0;
		}

		int count = 0;
		for (HgtEvent event : _events) {
			if (event.isBad()) {
				count++;
			}
		}

		return count;
	}

	/**
	 * This function counts the number of violated events. An event x -> y is
	 * violated if x is an ancestor, but not parent, of y.
	 *
	 * @return: The number of violated events.
	 */
	public int countViolatedEvents()
	{
		if (_events == null) {
			return 0;
		}

		int count = 0;
		for (HgtEvent event : _events) {
			if (event.isViolated()) {
				count++;
			}
		}

		return count;
	}

	/**
	 * This function returns a list of events within this scenario.
	 *
	 * @return: A list of events.
	 */
	public List<HgtEvent> getEvents()
	{
		return _events;
	}

	/**
	 * This function adds a new event to the list of events contained in this scenario.
	 *
	 * @param e: The event we want to add to this scenario.
	 *
	 * @return: <code>true</code> if this event is successfully added (i.e., it can be added
	 * to the list _events without resulting in any duplication); <code>false</code> otherwise.
	 */
	public boolean addEvent(HgtEvent e)
	{
		if (_events.contains(e)) {
			return false;
		}
		else {
			_events.add(e);
			return true;
		}
	}

	/**
	 * Add all events from another scenario to this scenario.
	 */
	public void addEvents(HgtScenario hs) {
		for (HgtEvent he : hs.getEvents()) {
			addEvent(he);
		}
	}

	/**
	 * This function removes an event from the list _events of this scenario.
	 *
	 * @param e: The event to be removed from this scenario.
	 *
	 * @return: <code>true</code> if the event is indeed contained in the scenario, and
	 * we can remove it; <code>false</code> otherwise.
	 */
	public boolean removeEvent(HgtEvent e)
	{
		int i = _events.indexOf(e);

		if (i == -1) {
			return false;
		}
		else {
			_events.remove(i);
			return true;
		}
	}

	/* Data members of ExHgtScenario */
	private List<HgtEvent> _events;	// List of HGT events; Each event is an object of HgtEvent.
}

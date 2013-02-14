package edu.rice.cs.bioinfo.programs.soranus.views.swing;

import y.view.Graph2DView;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/7/13
 * Time: 2:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class ZoomableGraph2DView extends JPanel
{
    private Graph2DView _view;

    ZoomableGraph2DView(Graph2DView view)
    {
        _view = view;
        this.setLayout(new BorderLayout());
        this.add(_view, BorderLayout.CENTER);

        JPanel toolbar = new JPanel();
        this.add(toolbar, BorderLayout.NORTH);

        JButton fitZoom = new JButton("Standard Zoom");
        fitZoom.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                _view.fitContent();
                _view.updateView();
            }
        });

        toolbar.add(fitZoom);

        JButton zoomIn = new JButton("Zoom In");
        zoomIn.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                _view.setZoom(_view.getZoom() * 1.3);
                _view.updateView();
            }
        });

        toolbar.add(zoomIn);

        JButton zoomOut = new JButton("Zoom Out");
        zoomOut.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                _view.setZoom(_view.getZoom() / 1.3);
                _view.updateView();
            }
        });

        toolbar.add(zoomOut);
    }
}

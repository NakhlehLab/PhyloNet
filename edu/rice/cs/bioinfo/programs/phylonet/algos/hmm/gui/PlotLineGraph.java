package edu.rice.cs.bioinfo.programs.phylonet.algos.hmm.gui;


import com.google.gson.Gson;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;

public class PlotLineGraph extends SavableTabData {

    JFreeChart chart;
    String title;
    String XLABEL;
    String YLABEL;
    double[] inputData;

    /**
     * Deserializes a PlotAreaGraph from a JsonObject
     *
     * @param obj The JsonObject holding the data.
     */
    public PlotLineGraph(JsonObject obj) {
        this(obj.get("title").getAsString(), obj.get("xlabel").getAsString(), obj.get("ylabel").getAsString(), new Gson().fromJson(obj.get("values"), double[].class));
    }

    /**
     * Creates a PlotAreaGraph from data with a given title.
     *
     * @param title     The title of the graph.
     * @param inputData The data, where inputData[x] is the y value for a type at a given x coordinate.
     */
    public PlotLineGraph(String title, String XLabel, String YLabel, double[] inputData) {
        this.title = title;
        this.XLABEL = XLabel;
        this.YLABEL = YLabel;
        this.inputData = inputData;

        XYDataset dataset = createDataset(inputData);
        chart = createChart(title, XLabel, YLabel, dataset);
    }


    /**
     * Creates a dataset from a two dimensional array.
     *
     * @param inputData The data, where inputData[x][type] is the y value for a type at a given x coordinate.
     * @return A dataset holding the given data.
     */
    private XYDataset createDataset(double[] inputData) {
        XYSeriesCollection collection = new XYSeriesCollection();

        collection.addSeries(createSeries("", inputData));
        return collection;
    }

    /**
     * Creates a series from data with a given label.
     *
     * @param label   The label for the series.
     * @param doubles The data for the series.
     * @return The corresponding XYSeries.
     */
    private XYSeries createSeries(Comparable label, double[] doubles) {
        XYSeries series = new XYSeries(label, false, false);
        for (int i = 0; i < doubles.length; i++) {
            series.add(i, doubles[i]);
        }
        return series;
    }


    /**
     * Creates a stacked XY area chart with the given title and dataset.
     *
     * @param title   The title of the graph.
     * @param dataset The dataset to display
     * @return The chart.
     */
    private JFreeChart createChart(String title, String XAxisLabel, String YAxisLabel, XYDataset dataset) {

        JFreeChart chart = ChartFactory.createScatterPlot(
                title, XAxisLabel, YAxisLabel, dataset, PlotOrientation.VERTICAL, true, true, false);

        chart.setBackgroundPaint(Color.white);

        final XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.lightGray);
        plot.setRangeGridlinePaint(Color.white);

        // customise the range axis...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setAutoRangeIncludesZero(false);

        final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
        renderer.setBaseShapesVisible(false);
        Stroke stroke = new BasicStroke(4);
        renderer.setSeriesStroke(0,stroke);
        plot.setRenderer(renderer);

        return chart;
    }

    @Override
    public JComponent getComponent() {
        return new ChartPanel(chart);
    }

    @Override
    public String getType() {
        return "PlotLineGraph";
    }

    @Override
    public JsonElement getData() {

        Gson gson = new Gson();

        JsonObject obj = new JsonObject();
        obj.addProperty("title", title);
        obj.addProperty("xlabel", XLABEL);
        obj.addProperty("ylabel", YLABEL);
        obj.add("values", gson.toJsonTree(inputData));

        return obj;
    }

}
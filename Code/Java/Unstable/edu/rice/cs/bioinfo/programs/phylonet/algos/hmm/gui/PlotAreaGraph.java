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
import org.jfree.data.xy.DefaultTableXYDataset;
import org.jfree.data.xy.TableXYDataset;
import org.jfree.data.xy.XYSeries;

import javax.swing.*;
import java.awt.*;

public class PlotAreaGraph extends SavableTabData {

    JFreeChart chart;
    String title;
    double[][] inputData;

    /**
     * Deserializes a PlotAreaGraph from a JsonObject
     *
     * @param obj The JsonObject holding the data.
     */
    public PlotAreaGraph(JsonObject obj) {
        this(obj.get("title").getAsString(), new Gson().fromJson(obj.get("values"), double[][].class));
    }

    /**
     * Creates a PlotAreaGraph from data with a given title.
     *
     * @param title     The title of the graph.
     * @param inputData The data, where inputData[x][type] is the y value for a type at a given x coordinate.
     */
    public PlotAreaGraph(String title, double[][] inputData) {
        this.title = title;
        this.inputData = inputData;

        TableXYDataset dataset = createDataset(inputData);
        chart = createChart(title, dataset);
    }


    /**
     * Creates a dataset from a two dimensional array.
     *
     * @param inputData The data, where inputData[x][type] is the y value for a type at a given x coordinate.
     * @return A dataset holding the given data.
     */
    private TableXYDataset createDataset(double[][] inputData) {
        DefaultTableXYDataset collection = new DefaultTableXYDataset();

        for (int j = 0; j < inputData.length; j++) {
            double[] doubles = inputData[j];
            XYSeries series = createSeries(j, doubles);
            collection.addSeries(series);
        }
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
    private JFreeChart createChart(String title, TableXYDataset dataset) {

        JFreeChart chart = ChartFactory.createStackedXYAreaChart(
                title, "Nucleotide", "Probability", dataset, PlotOrientation.VERTICAL, true, true, false);

        chart.setBackgroundPaint(Color.white);

        final XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.lightGray);
        plot.setRangeGridlinePaint(Color.white);

        // customise the range axis...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        rangeAxis.setAutoRangeIncludesZero(true);

        final NumberAxis domainAxis = (NumberAxis) plot.getDomainAxis();
        domainAxis.setRange(0,dataset.getItemCount());

        return chart;
    }

    @Override
    public JComponent getComponent() {
        return new ChartPanel(chart);
    }

    @Override
    public String getType() {
        return "PlotAreaGraph";
    }

    @Override
    public JsonElement getData() {

        Gson gson = new Gson();

        JsonObject obj = new JsonObject();
        obj.addProperty("title", title);
        obj.add("values", gson.toJsonTree(inputData));

        return obj;
    }

}







package runHmm;

public class MemoryReport {
    private static final long MiB = 1024 * 1024;

    public static void report() {
        Runtime runtime = Runtime.getRuntime();

        System.err.println();

        long total = runtime.totalMemory() / MiB;
        long free = runtime.freeMemory() / MiB;
        long used = total - free;
        System.err.println("Total: " + total + " MiB; Free: " + free + " MiB; --> Used: " + used + " MiB.");

        long max = runtime.maxMemory() / MiB;
        System.err.println("Max: " + max + " MiB.");
        System.err.println();
    }
}

package containers;

public interface IArray<Type> {
    public Type get(int index);
    public void set(int index, Type type);
    public int getLength();
}

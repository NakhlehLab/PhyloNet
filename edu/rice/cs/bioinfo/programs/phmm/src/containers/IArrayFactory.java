package containers;

public interface IArrayFactory<Type> {
    public IArray<Type> make(int length);
}

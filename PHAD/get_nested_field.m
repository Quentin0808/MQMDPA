function val = get_nested_field(S, field1, field2, default_val)

    if isfield(S, field1) && isfield(S.(field1), field2)
        val = S.(field1).(field2);
    else
        val = default_val;
    end

end

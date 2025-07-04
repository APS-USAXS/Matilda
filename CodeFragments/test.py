def flatten_dict(d):
    items = []
    for k, v in d.items():
        if isinstance(v, dict):
            items.extend(flatten_dict(v).items())
        elif isinstance(v, bytes):  # Convert bytes to string
            items.append((k, v.decode('utf-8')))
        elif not isinstance(v, list):  # Skip arrays (lists)
            items.append((k, v))
    return dict(items)

def dict_to_string(d):
    flat_dict = flatten_dict(d)
    return ';'.join(f'{key}={value}' for key, value in flat_dict.items())

# Example usage:
nested_dict = {
    'monochromator': {
        'energy': '0.5th',
        'wavelength': 42,
        'inner_key3': [1, 2, 3]  # This will be skipped
    }
}


print(dict_to_string(nested_dict))


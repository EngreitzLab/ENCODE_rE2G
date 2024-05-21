import pickle
import click
import json


def get_params(default_params, override_params):
    final_params = default_params

    # replace overriding parameters
    if isinstance(override_params, dict):
        for arg, val in override_params.items():
            final_params[arg] = val

    # replace null with None
    final_params = {
        key: None if value == "null" else value for key, value in final_params.items()
    }

    # convert numeric strings to correct type
    for key, val in final_params.items():
        if isinstance(val, str):
            try:
                if key == "max_iter" or key == "random_state" or key == "n_jobs":
                    final_params[key] = int(float(val))
                else:
                    final_params[key] = float(val)
            except ValueError:
                pass
    return final_params


@click.command()
@click.option("--default_params", required=True)
@click.option("--override_params", required=True)
@click.option("--output_file", required=True)
def main(default_params, override_params, output_file):
    # read in default and override params as strings and convert to dict
    default_params_fixed = (
        default_params.replace("'", '"')
        .replace("True", "true")
        .replace("False", "false")
        .replace("None", "null")
    )
    override_params_fixed = (
        override_params.replace("'", '"')
        .replace("True", "true")
        .replace("False", "false")
        .replace("None", "null")
    )

    default_params_dict = json.loads(default_params_fixed)
    override_params_dict = json.loads(override_params_fixed)

    # generate final params
    final_params = get_params(default_params_dict, override_params_dict)

    # save to pickle
    with open(output_file, "wb") as f:
        pickle.dump(final_params, f)


if __name__ == "__main__":
    main()

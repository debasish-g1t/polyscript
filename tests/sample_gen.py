from polyscript.polymerizer import MpiPolymerizer


polymerizer = MpiPolymerizer(
    exp_name="sample_outputs",
    input_data_split_dir = "./test_data/cls_split_10",
    output_path="./test_data/",
)

polymerizer.run()

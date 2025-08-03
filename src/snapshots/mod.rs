use crate::spin::SpinState;
use ndarray::Array4;

pub fn save_snapshots_to_hdf5<S: SpinState + hdf5_metno::H5Type>(
    filename: &str,
    equil_data: &[Array4<S>],
    steps_data: &[Array4<S>],
) -> hdf5_metno::Result<()> {
    let file = hdf5_metno::File::create(filename)?;
    let equil_views: Vec<_> = equil_data.iter().map(|a| a.view()).collect();
    let equil_stacked = ndarray::stack(ndarray::Axis(0), &equil_views)?;
    let steps_views: Vec<_> = steps_data.iter().map(|a| a.view()).collect();
    let steps_stacked = ndarray::stack(ndarray::Axis(0), &steps_views)?;

    let _equil_ds = file
        .new_dataset_builder()
        .with_data(&equil_stacked)
        .deflate(9)
        .create("snapshots/equil")
        .unwrap();
    let _steps_ds = file
        .new_dataset_builder()
        .with_data(&steps_stacked)
        .deflate(9)
        .create("snapshots/steps")
        .unwrap();
    Ok(())
}

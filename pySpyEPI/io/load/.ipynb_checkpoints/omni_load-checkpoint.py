"""Load OMNI data
Created 19 February 2026
alanahco

List Functions
--------------
get_OMNI
"""
from cdasws import CdasWs, TimeInterval
from cdasws.datarepresentation import DataRepresentation as dr


def get_OMNI(t1, t2):
    """ obtain OMNI data and return as a dataframe for OMNI_HRO_5MIN

    Parameters
    ----------
    t1 : datetime
        starting time
    t2 : datetime
        ending time

    Returns
    -------
    OMNI_data : dataframe
        dataframe of OMNI 5 min data
    """
    cdas = CdasWs()
    datasets = cdas.get_datasets(observatoryGroup='ACE',
                                 instrumentType='Magnetic Fields (space)')
    for index, dataset in enumerate(datasets):
        dataset_id = dataset["Id"]

        if 'OMNI_HRO_5MIN' in dataset_id:  # select this dataset for use below
            break
    doi = dataset['Doi']

    # get times in proper string format
    t1s = t1.strftime('%Y-%m-%dT%H:%M:%SZ')
    t2s = t2.strftime('%Y-%m-%dT%H:%M:%SZ')
    variables = cdas.get_variables(doi)
    var_names = []
    for index, variable in enumerate(variables):
        name = variable['Name']
        var_names.append(name)

    time_interval = TimeInterval(t1s, t2s)

    # GET OMNI data for specified time interval
    _, data = cdas.get_data(doi, var_names, time_interval,
                            dataRepresentation=dr.XARRAY)
    OMNI_data = data[var_names].to_dataframe()

    return OMNI_data

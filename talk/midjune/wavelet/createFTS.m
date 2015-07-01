function fts = createFTS(T)
    dates = datenum(table2array(T(:,1)),'dd/mm/yyyy');
    data = T.Value;
    fts = fints(dates, data);
end
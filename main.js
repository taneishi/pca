$(function(){
    var style = 'border: 1.5px solid #dfdfdf; padding: 10px;';

    $('#container').w2layout({
        name: 'container',
        panels: [
            { type: 'left', size: 700, style: style, resizable: true },
            { type: 'main', style: style, content: '<canvas id="pca" width="600px" height="600px"></canvas>' },
        ]
    });

    d3.csv('https://raw.githubusercontent.com/taneishi/pca/master/data/iris.csv').then(function(data){
        // rename column names
        X = data.map(row => data.columns.slice(0, 4).map(col => parseFloat(row[col])));
        y = data.map(row => row['Species']);

        pc = new PCA().pca(X, 2); // R-compatible

        labels = Array.from(new Set(y));
        colors = ['#aea', '#aae', '#eaa']

        datasets = labels.map(function(label, i){
            return {
                label: label, 
                data: pc.filter((row, j) => (y[j] == label)),
			    borderColor: '#222',
                backgroundColor: colors[i],
                pointRadius: 5,
                };
        });

        columns = data.columns.map(col => col.replace('.', ' '));
        columns = columns.concat(['PC1', 'PC2']);

        data = X.map(function(row, i){
            row = row.concat(y[i]).concat(pc[i]);
            row = row.map((val, j) => [columns[j], val]);
            row = [['recid',i]].concat(row)
            return Object.fromEntries(row);
        });

        columns = columns.map(col => ({ field: col, text: col, size: '25%', sortable: true, resizable: true }));

        $('#layout_container_panel_left').w2grid({
            name: 'grid',
            header: 'Iris dataset',
            show: {
                toolbar: true,
                footer: true,
                lineNumbers: true,
            },
            columns: columns,
            searches: [
                { field: 'Species', text: 'Species', type: 'text' },
            ],
            sortData: [{ field: 'Species', direction: 'ASC' }],
            records: data,
        });

        var mychart = new Chart('pca', {
            type: 'scatter',
            data: { datasets: datasets },
            options: {
                responsive: false,
                scales: {
                    x: {
                        title: { display: true, text: 'PC1', },
                    },
                    y: {
                        title: { display: true, text: 'PC2', },
                    },
                },
                plugins: {
                    legend: {
                        position: 'top',
                    },
                    title: {
                        display: true,
                        text: 'Principal Component Anaylsis'
                    },
                }
            },
        });

    })
});

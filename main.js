$(function(){
    var style = 'border: 1.5px solid #dfdfdf; padding: 10px;';

    $('#container').w2layout({
        name: 'container',
        panels: [
            { type: 'left', size: 720, style: style, resizable: true },
            { type: 'main', style: style, content: '<canvas id="pca" width="600px" height="600px"></canvas>' },
        ]
    });


    d3.csv('https://raw.githubusercontent.com/taneishi/pca/master/data/iris.csv').then(function(data){
        columns = data.columns;

        // rename column names
        data = data.map(function(row){
            columns.slice(0, 4).map(function(col){
                row = Object.assign(row, {[col.replace('.', ' ')]: Number(row[col])});
                delete row[col];
            });
            return row;
        });

        columns = columns.map(col => col.replace('.', ' '));

        X = data.map(row => 
            columns.slice(0, 4).map(col => row[col])
        );

        pc = new PCA().pca(X, 2);

        labels = data.map(row => row['Species']);
        pc = pc.map((row, i) => row.concat(labels[i]));

        labels = Array.from(new Set(labels));

        colors = ['#aea', '#aae', '#eaa']

        datasets = labels.map(function(label){
            return {
                label: label, 
                data: pc.filter(row => (row[4] == label)),
			    borderColor: '#222',
                backgroundColor: colors[labels.indexOf(label)],
                pointRadius: 5,
                };
        });

        columns = columns.map(col => ({ field: col, text: col, size: '25%', sortable: true, resizable: true }));
        data = data.map((row, i) => Object.assign({recid: i}, row));

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

